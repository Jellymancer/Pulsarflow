#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool } from './modules'
include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'
include { generateDMFiles as generateDMFiles } from './modules'
include { calculate_max_acceleration as calculate_max_acceleration } from './modules'
include { peasoup as peasoup } from './modules'
include {split_filterbank as split_filterbank} from './modules'
include {fold_peasoup_cands_pulsarx as fold_peasoup_cands_pulsarx} from './modules'
include {aggregate_peasoup_output as aggregate_peasoup_output} from './modules'
include {classify_candidates as classify_candidates} from './modules'
include {create_tar_archive as create_tar_archive} from './modules'
include {sift_candidates as sift_candidates} from './modules'



workflow {

    filterbank_channel_with_metadata = Channel
        .fromPath("${params.filterbank_list}")
        .splitCsv(header: true, sep:',') // Ensure the separator is specified if not comma
        .map { row ->
            // Now each row is a Groovy map with column names as keys
            def filterbank_files = row.filterbank_files.trim() // Trim leading and trailing spaces
            def target = row.target.trim()
            def beam_id = row.beam_num.trim()
            def pointing_id = row.pointing_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-") // Replace space with dash in utc_start
            return tuple(filterbank_files, target, pointing_id, beam_id, utc_start)
        }

    // Process the filterbank file with filtool
    if (params.use_filtool == 1) {
        processed_filterbank = filtool(filterbank_channel_with_metadata, params.filtool_rfi_filter_list,
        params.filtool_rfi_channel_mask_pairs, params.filtool_threads, params.telescope)
    }

    // Update the filterbank channel in case filtool is used
    if (params.use_filtool != 1) {
        updated_filterbank_channel = filterbank_channel_with_metadata
    } else {
        updated_filterbank_channel = processed_filterbank.map { metadata, filepath ->
            def (raw_filterbank, target, pointing_id, beam_name, utc_start) = metadata
            // Trim the string values
            target = target.trim()
            pointing_id = pointing_id.trim()
            beam_name = beam_name.trim()
            utc_start = utc_start.trim()
            return tuple(filepath, target, pointing_id, beam_name, utc_start)
        }
    }

    // Split the filterbank files if enabled
    if (params.split_filterbank == 1) {
        split_channel = split_filterbank(updated_filterbank_channel, params.split_blocks ?: 2)
        
        // Split the filterbank output into individual filterbank channels
        split_filterbank_files = split_channel
            .flatMap { file, target_name, pointing_id, beam_name, utc ->
                file.collect { f ->
                    // Get the base filename (without extension)
                    def fileName = f.getBaseName()

                    // Split the filename by a delimiter (e.g., '_') and take the first part
                    def split_id = fileName.split('_')[0]

                    // Return a tuple with the file and generated ID
                    return tuple(f, split_id, target_name, pointing_id, beam_name, utc)
                }
            }
    } else
        // If splitting is not enabled, just use the updated filterbank channel, setting the split_id to 0
        split_filterbank_files = updated_filterbank_channel.map {filepath, target, pointing_id, beam_name, utc_start ->
            return tuple(filepath, 0, target, pointing_id, beam_name, utc_start)
        }

    
    // Generate DM files
    dm_files = generateDMFiles()
    dm_files_channel = Channel.fromPath("${params.dm_out_dir}/*.txt")
    .map { dm_file ->
        def dm_file_name = dm_file.getBaseName()
        return tuple(dm_file, dm_file_name)
    }

    // Find the nearest power of two
    nearest_two_output = nearest_power_of_two_calculator_apsuse(split_filterbank_files)

    // Set an acceleration limit
    if (params.calculate_max_accel == 1){
        // Automatically calculate the max acceleration range based on the filterbank length
        max_acceleration_output = calculate_max_acceleration(nearest_two_output)
    }
    else {
        // If the max acceleration is manually set, use the provided value
        max_acceleration_output = nearest_two_output.map { target_name, pointing_id, split_id, beam_name, utc, nearest_power_of_2 ->
            return tuple(target_name, pointing_id, split_id, beam_name, utc, nearest_power_of_2, params.acc_start, params.acc_end)
        }
    }
    // max_acceleration_output.view()  
    // Launch peasoup
    peasoup_channel = max_acceleration_output.combine(dm_files_channel)
    // peasoup_channel.view()
    peasoup_output = peasoup(peasoup_channel, params.total_cands_limit, params.min_snr, params.ram_limit_gb, params.nh, params.ngpus, params.kill_file)

    unique_fil_name = peasoup_output
        .map { target_name, pointing_id, split_id, dm_range, beam_name, utc, dm_file, fil_file, xml_file ->
            // Extract the base name of the fil_file
            def fil_base_name = fil_file.getBaseName()

            // Return a tuple with the base name and the original file path for grouping
            return tuple(fil_base_name, fil_file)
        }
        // Group by the 5th element of the tuple, which is fil_base_name
        .groupTuple(by: 0)
        .map { fil_base_name, fil_files->
            // fil_files is a list now, we select the first one
            def first_fil_file = fil_files[0]

            // Return a flat tuple without the fil_files list
            return tuple(fil_base_name, first_fil_file)
        }

    peasoup_no_filname = peasoup_output
        .map { target_name, pointing_id, split_id, dm_range, beam_name, utc, dm_file, fil_file, xml_file  ->

            def fil_base_name = fil_file.getBaseName()
            return tuple(fil_base_name, target_name, pointing_id, split_id, dm_range, beam_name, utc, dm_file, fil_file, xml_file)
        }
    //
    // Aggreagate the peasoup output based on pointing_id and beam name. 
    // This essentially aggregates all candidates from the different DM subdivisions
    grouped_peasoup_output = peasoup_no_filname.combine(unique_fil_name, by: 0)
    .map {fil_base_name, target_name, pointing_id, split_id, dm_range, beam_name, utc, dm_file, fil_file, xml_file, first_fil_file->
        
        return tuple(target_name, pointing_id, split_id, dm_range, beam_name, utc, first_fil_file, xml_file)
    }.groupTuple(by: [0, 1, 2, 5])
    // Convert the xmls into a single cvs file
    aggregated_peasoup_output = aggregate_peasoup_output(grouped_peasoup_output)

    // Sifting
    sifted_candidates = sift_candidates(aggregated_peasoup_output)
    // sifted_candidates.view()
    to_fold = aggregated_peasoup_output
    .flatMap { tuple ->
        def(target_name, pointing_id, split_id, dm_range, beam_names, utc, fil_files, cand_file, meta_file) = tuple

        // Convert `fil_files` and `dm_range` to lists if they are not already lists
        fil_files = fil_files instanceof List ? fil_files : [fil_files]
        dm_range = dm_range instanceof List ? dm_range : [dm_range]
        beam_names = beam_names instanceof List ? beam_names : [beam_names]

        def result = []
        for (int i = 0; i < fil_files.size(); i++) {
            result << [target_name, pointing_id, split_id, dm_range[i], beam_names[i], utc, fil_files[i], cand_file, meta_file ]
        }
        // fil_files.zip(beam_names).collect { fil_file, beam_name -> [target_name, pointing_id, split_id, beam_name, utc, fil_file, cand_file, meta_file] }
    return result
    }

    // Pulsarx folding
    pulsarx_output = fold_peasoup_cands_pulsarx(to_fold)

    // Now, we group all of the
    grouped_by_beam = pulsarx_output.groupTuple(by: [0, 1, 2])
    .map{target_name, pointing_id, split_id, dm_range, beam_names, utc, meta_files, archive_files, png_files, candfiles, cands_files ->
        def flat_archive_files = archive_files.flatten()
        def flat_png_files = png_files.flatten()
        def flat_candfiles = candfiles.flatten()
        def flat_cands_files = cands_files.flatten()
        def meta_file = meta_files[0]

        return tuple(target_name, pointing_id, split_id, dm_range, beam_names, utc, meta_file, flat_archive_files, flat_png_files, flat_candfiles, flat_cands_files)
    }
    grouped_by_beam.view()

    classified_candidates = classify_candidates(grouped_by_beam)

    tar_output = create_tar_archive(classified_candidates, params.metafile)
    

}