#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool } from './modules'

include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'

include { generateDMFiles as generateDMFiles } from './modules'


process split_filterbank {
    label 'split_filterbank'
    container "${params.presto5_singularity_image}"
    publishDir "out/SPLIT/${utc}/${target_name}/${beam_name}/", pattern: "*.fil", mode: 'symlink'

    input:
    tuple path(fil_file), val(target_name), val(pointing_id), val(beam_name), val(utc)
    val split_blocks

    output:
    tuple path("*.fil"), val(target_name), val(pointing_id), val(beam_name), val(utc)

    script:
    """
    python3 "${params.split_script}" ${fil_file} --nblocks ${split_blocks}
    """
}

process peasoup {
    label 'peasoup'
    container "${params.search_singularity_image}"

    publishDir { "out/SEARCH/${utc}/${target_name}/${beam_name}/${pointing_id}/${dm_file.baseName}" }, pattern: "*.xml", mode: 'copy'
    input:
    tuple path(fil_file), val(split_id), val(target_name), val(pointing_id), val(beam_name), val(utc), val(fft_size), path(dm_file), val(dm_range)
    val(total_cands_limit)
    val(min_snr)
    val(acc_start)
    val(acc_end)
    val(ram_limit_gb)
    val(nh)
    val(ngpus)
    val(kill_file)

    output:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(dm_file), path(fil_file, followLinks: false), path("*.xml")

    script:
    """
    #!/bin/bash

    if [ -z "$kill_file" ]; then
        kill_file_option=\${kill_file ? "-k \${kill_file}" : ""}
    else
        kill_file_option=""
    fi

    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus} \$kill_file_option
    
    # Rename the output file
    mv **/*.xml ${beam_name}_${dm_file.baseName}_overview.xml
    """
}


process aggregate_peasoup_output {
    label 'aggregate_peasoup_output'
    container "${params.utility_singularity_image}"
    publishDir "out/SEARCH/${utc}/${target_name}/${split_id}/${pointing_id}", pattern: "*.{csv, meta}", mode: 'copy'

    stageInMode 'symlink'

    input:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_file, stageAs: "?/*"), path(xml_files)

    output:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_file),  path("*.csv"), path("*.meta")

    script:
    // join input on pointing_id
    """
    echo ""
    python3 "${params.aggregate_script}"  ${xml_files}
    """
}

process sift_candidates {
    container "${params.utility_singularity_image}"
    publishDir "out/SIFTING/${utc}/${target_name}/${split_id}/${pointing_id}/", pattern: "*.csv", mode: 'symlink'

    input:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_files, stageAs: "?/*"), path(cand_file), path(meta_file)

    output:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_files), path("*.csv"), path(meta_file)

    script:
    """
    python3 ${params.sift_script} -c  /hercules/results/jjawor/GC/NGC2808/subband_followup/cfbf_test/sifting/default_config.json -m ${meta_file} ${cand_file}
    """
}

process fold_peasoup_cands_pulsarx {
    label 'pulsarx'
    container "${params.pulsarx_singularity_image}"
    publishDir "out/FOLDING/${utc}/${target_name}/${split_id}/${pointing_id}/${beam_name}/", pattern: "*.{ar,png,xml,candfile,cands}", mode: 'symlink'

    input:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_file), path(cand_file), path(meta_file)

    output:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(meta_file), path("*.ar"), path("*.png"), path("*.candfile"), path("*.cands")

    script:
    """
    python3 ${params.fold_script} -i ${cand_file} -if ${fil_file} -t pulsarx -p ${params.pulsarx_fold_template} -b ${beam_name} -threads ${params.psrfold_fil_threads} -ncands ${params.no_cands_to_fold} -n ${params.psrfold_fil_nh} --snr_min ${params.min_snr_fold}  --metafile ${meta_file}
    
    # Get the base name of the fil_file
    echo "qwer"
    CANDSFILE=*.cands
    BASECANDSFILE=\$(basename \${CANDSFILE})
    NEWCANDFILE="\${BASECANDSFILE%.*}".candfile
    mv pulsarx.candfile \${NEWCANDFILE}

    """

}

process classify_candidates {
    label 'classify_candidates'
    container "${params.old_pulsarx_singularity_image}"
    publishDir "out/CLASSIFICATION/${target_name}/${split_id}/${pointing_id}/", pattern: "*.csv", mode: 'symlink'

    input:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(meta_file), path(archives), path(pngs), path(candfiles), path(cands)

    output:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(archives), path(pngs), path(candfiles), path(cands), path("*.csv")

    script:
    """
    python2 "${params.classify_script}" -m ${params.model_dir}
    """
}

process create_tar_archive {
    label 'create_tar_archive'
    container "${params.pulsarx_singularity_image}"
    publishDir "out/TAR/${target_name}/${split_id}/${pointing_id}/", pattern: "out/*.tar", mode: 'copy'

    input:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(archives), path(pngs), path(candfiles), path(cands), path(full_candfile)
    path(metafile)

    output:
    path("out/*.tar")

    script:
    """
    python3 ${params.prepare_for_candyjar_script} -d \${PWD} -m ${metafile} -p ${pointing_id} -o \${PWD}/out 
    """
}

workflow {

    filterbank_channel_with_metadata = Channel
        .fromPath("${params.filterbank_list}")
        .splitCsv(header: true, sep:',') // Ensure the separator is specified if not comma
        .map { row ->
            // Now each row is a Groovy map with column names as keys
            def filterbank_files = row.filterbank_files.trim() // Trim leading and trailing spaces
            def target = row.target.trim()
            def beam_num = row.beam_num.trim()
            def pointing_id = row.pointing_id.trim()
            def utc_start = row.utc_start.trim().replace(" ", "-") // Replace space with dash in utc_start
            return tuple(filterbank_files, target, pointing_id, beam_num, utc_start)
        }

    // Process the filterbank file with filtool
    if (params.use_filtool == 1) {
        processed_filterbank = filtool(filterbank_channel_with_metadata, params.filtool_rfi_filter, params.filtool_threads, params.telescope, params.filtool_channel_mask)
    }

    // Check if neither processing flags are set to 1
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

    split_filterbank_files.view()
    // Find the nearest power of two
    nearest_two_output = nearest_power_of_two_calculator_apsuse(split_filterbank_files)
    
    // Launch peasoup
    peasoup_channel = nearest_two_output.combine(dm_files_channel)
    peasoup_output = peasoup(peasoup_channel, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.kill_file)

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

    grouped_peasoup_output = peasoup_no_filname.combine(unique_fil_name, by: 0)
    .map {fil_base_name, target_name, pointing_id, split_id, dm_range, beam_name, utc, dm_file, fil_file, xml_file, first_fil_file->
        
        return tuple(target_name, pointing_id, split_id, dm_range, beam_name, utc, first_fil_file, xml_file)
    }.groupTuple(by: [0, 1, 2, 5])

aggregated_peasoup_output = aggregate_peasoup_output(grouped_peasoup_output)
// aggregated_peasoup_output.view()

    // Aggregate based Sift the candidates
    // aggregated_peasoup_output.collect().set { collected_peasoup_output }
    // sifted_candidates = sift_candidates(aggregated_peasoup_output)
    // aggregated_peasoup_output.view()
    to_fold = aggregated_peasoup_output
    .flatMap { tuple ->
        def(target_name, pointing_id, split_id, dm_range, beam_names, utc, fil_files, cand_file, meta_file) = tuple

        def result = []
        for (int i = 0; i < fil_files.size(); i++) {
            result << [target_name, pointing_id, split_id, dm_range[i], beam_names[i], utc, fil_files[i], cand_file, meta_file ]
        }
        // fil_files.zip(beam_names).collect { fil_file, beam_name -> [target_name, pointing_id, split_id, beam_name, utc, fil_file, cand_file, meta_file] }
    return result
    }
    // // Pulsarx folding
    pulsarx_output = fold_peasoup_cands_pulsarx(to_fold)
//     // Classify candidates and create a TAR archive
// val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path("*.ar"), path("*.png"), path("*.candfile"), path("*.cands")
    grouped_by_beam = pulsarx_output.groupTuple(by: [0, 1, 2])
    .map{target_name, pointing_id, split_id, dm_range, beam_names, utc, meta_files, archive_files, png_files, candfiles, cands_files ->
        def flat_archive_files = archive_files.flatten()
        def flat_png_files = png_files.flatten()
        def flat_candfiles = candfiles.flatten()
        def flat_cands_files = cands_files.flatten()
        def meta_file = meta_files[0]

        return tuple(target_name, pointing_id, split_id, dm_range, beam_names, utc, meta_file, flat_archive_files, flat_png_files, flat_candfiles, flat_cands_files)
    }
// grouped_by_beam.view()

    classified_candidates = classify_candidates(grouped_by_beam)

    tar_output = create_tar_archive(classified_candidates, params.metafile)
    

}