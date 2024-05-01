#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool } from './modules'

include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'

include { generateDMFiles as generateDMFiles } from './modules'


process peasoup {
    label 'peasoup'
    container "${params.search_singularity_image}"

    publishDir { "/SEARCH/${dm_file.baseName}/" }, pattern: "**/*.xml", mode: 'copy'

    input:
    path(dm_file) 
    path(fil_file)
    val(target_name)
    val(beam_name)
    val(utc)
    val(fft_size)
    val(total_cands_limit)
    val(min_snr)
    val(acc_start)
    val(acc_end)
    val(ram_limit_gb)
    val(nh)
    val(ngpus)
    val(kill_file)

    output:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size), path("**/*.xml")

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
    // publishDir "RESULTS/${POINTING}/${UTC_OBS}/${BAND}/${BEAM}/04_FOLDING/", pattern: "*.ar", mode: 'copy'
    // publishDir "RESULTS/${POINTING}/${UTC_OBS}/${BAND}/${BEAM}/04_FOLDING/", pattern: "*.png", mode: 'copy'

    input:
    tuple path(input_file), val(POINTING), val(BEAM), val(UTC_OBS), val(fft_size), path(peasoup_xml_out)
    val(psrfold_fil_threads) 
    val(no_cands_to_fold)
    val(cmask)

    output:
    path("*.ar")
    path("*.png")

    script:
    """
    python3 ${params.fold_script} -i ${peasoup_xml_out} -t pulsarx -p ${params.pulsarx_fold_template} -b ${BEAM} -threads ${psrfold_fil_threads} -ncands ${no_cands_to_fold} -c ${cmask}
    """

}


workflow {
    // Process the filterbank file with filtool_apsuse
    if (params.use_filtool == 1){
        processed_filterbank = filtool(params.filterbank_file, params.target_name, params.beam_name, params.utc_start, params.filtool_rfi_filter, params.filtool_threads, params.telescope, params.filtool_channel_mask)
    }
    // Generate DM files
    nearest_two_output = nearest_power_of_two_calculator_apsuse(processed_filterbank)
    dm_file_path = generateDMFiles(params.output_dm_dir)
    all_dm_files = Channel.fromPath("${dm_file_path}/*.txt")
    peasoup_output = peasoup(all_dm_files, processed_filterbank, params.target_name, params.beam_name, params.utc_start, nearest_two_output, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.kill_file)

    // Pulsarx processing
    pulsarx_output = fold_peasoup_cands_pulsarx(peasoup_output, params.psrfold_fil_threads, params.no_cands_to_fold, params.psrfold_fil_channel_mask)

}
