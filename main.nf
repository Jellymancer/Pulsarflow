#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool } from './modules'

include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'

include { generateDMFiles as generateDMFiles } from './modules'





process peasoup {
    label 'peasoup'
    container "${params.search_singularity_image}"

    publishDir { "/SEARCH/${utc}/${target_name}/${beam_name}/${dm_file.baseName}" }, pattern: "**/*.xml", mode: 'copy'

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size), path(dm_file) 
    val(total_cands_limit)
    val(min_snr)
    val(acc_start)
    val(acc_end)
    val(ram_limit_gb)
    val(nh)
    val(ngpus)
    val(kill_file)

    output:
    tuple val(beam_name), val(utc), path("**/*.xml")

    script:
    """
    #!/bin/bash

    if [ -z "$kill_file" ]; then
        kill_file_option=\${kill_file ? "-k \${kill_file}" : ""}
    else
        kill_file_option=""
    fi

    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus} \$kill_file_option
    """
}


process fold_peasoup_cands_pulsarx {
    label 'pulsarx'
    container "${params.pulsarx_singularity_image}"
    publishDir "FOLDING/${utc}/${target_name}/${beam_name}/${dm_file.baseName}/", pattern: "*.{ar,png,xml,candfile,cands}", mode: 'symlink'

    input:
    tuple val(BEAM), val(UTC_OBS), path(xml_file)

    output:
    tuple path("*.ar"), path("*.png")

    script:
    """
    python3 ${params.fold_script} -i ${xml_file} -t pulsarx -p ${params.pulsarx_fold_template} -b ${BEAM} -threads ${params.psrfold_fil_threads} -ncands ${params.no_cands_to_fold} -c ${params.cmask}
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
            def utc_start = row.utc_start.trim().replace(" ", "-") // Replace space with dash in utc_start
            return tuple(filterbank_files, target, beam_num, utc_start)
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
            def (raw_filterbank, target, beam_name, utc_start) = metadata
            // Trim the string values
            target = target.trim()
            beam_name = beam_name.trim()
            utc_start = utc_start.trim()
            return tuple(filepath, target, beam_name, utc_start)
        }
    }

    // Generate DM files
    dm_files = generateDMFiles()
    dm_files_channel = Channel.fromPath("${params.dm_out_dir}/*.txt")
    // dm_files_channel.view()

    // Find the nearest power of two
    nearest_two_output = nearest_power_of_two_calculator_apsuse(updated_filterbank_channel)

    // Launch peasoup
    peasoup_channel = nearest_two_output.combine(dm_files_channel)
    peasoup_output = peasoup(peasoup_channel, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.kill_file)

    // Pulsarx folding
    pulsarx_output = fold_peasoup_cands_pulsarx(peasoup_output)

}
