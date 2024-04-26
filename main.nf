#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//In nextflow dsl2 syntax, you cannot re-use the same process multiple times. Instead we make modules and call them separate names for APSUSE and PTUSE to solve the corner case where you search and fold both PTUSE & APSUSE observations.

include { filtool as filtool_apsuse } from './modules'
include { filtool_ptuse } from './modules'

include { digifil as digifil_apsuse } from './modules'
include { digifil_ptuse } from './modules'



include { create_phase_predictor as create_phase_predictor_apsuse } from './modules'
include { create_phase_predictor as create_phase_predictor_ptuse } from './modules'

include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_apsuse } from './modules'
include { nearest_power_of_two_calculator as nearest_power_of_two_calculator_ptuse } from './modules'


//include { dspsr_fold_phase_predictor as apsuse_fold_phase_predictor } from './modules'
include { dspsr_fold_phase_predictor_parallel as apsuse_fold_phase_predictor_parallel } from './modules'
include { dspsr_fold_phase_predictor_serial as apsuse_fold_phase_predictor_serial } from './modules'

include { dspsr_fold_phase_predictor as ptuse_fold_phase_predictor } from './modules'

include { dspsr_fold_ephemeris as apsuse_fold_ephemeris } from './modules'
//include { dspsr_fold_ephemeris as ptuse_fold_ephemeris } from './modules'
include { dspsr_fold_ephemeris_ptuse as ptuse_fold_ephemeris } from './modules'

include { clfd as clfd_apsuse_predictor } from './modules'
include { clfd as clfd_apsuse_eph } from './modules'
include { clfd as clfd_ptuse_predictor } from './modules'
include { clfd as clfd_ptuse_eph } from './modules'

include { pdmp as pdmp_apsuse_predictor } from './modules'
include { pdmp as pdmp_apsuse_eph } from './modules'
include { pdmp as pdmp_ptuse_predictor } from './modules'
include { pdmp as pdmp_ptuse_eph } from './modules'



process query_db {
    label 'query_db'
    container "${params.sql_query_image}"
    input:
    val(target_name)
    val(beam_name)
    val(utc_start)
    val(utc_end)

    output:
    path "grouped_files.txt"

    script:
    def beamArg = beam_name ? "-b ${beam_name}" : ""
    def utc_arg = utc_start ? "-us ${utc_start} -ue ${utc_end}" : ""
    """
    python ${params.db_query} -t ${target_name} ${beamArg} ${utc_arg}
    """

}

process peasoup_apsuse {
    label 'peasoup'
    container "${params.search_singularity_image}"


    publishDir "SEARCH/${params.target_name}/APSUSE/${utc}/${beam_name}/${params.dm_file.split('\\.')[0]}/", pattern: "**/*.xml", mode: 'copy'

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size)
    val(dm_file) 
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
   
    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus} -k ${kill_file}
    """
}


process peasoup_ptuse {
    label 'peasoup'
    container "${params.search_singularity_image}"
    publishDir "SEARCH/${params.target_name}/PTUSE/${utc}/${beam_name}/", pattern: "**/*.xml", mode: 'copy'
    

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc), val(fft_size)
    path(dm_file) 
    val(total_cands_limit)
    val(min_snr)
    val(acc_start)
    val(acc_end)
    val(ram_limit_gb)
    val(nh)
    val(ngpus)
    

    output:
    path("**/*.xml")

    script:
    """
    #!/bin/bash
   
    peasoup -i ${fil_file} --fft_size ${fft_size} --limit ${total_cands_limit} -m ${min_snr} --acc_start ${acc_start} --acc_end ${acc_end} --dm_file ${dm_file} --ram_limit_gb ${ram_limit_gb} -n ${nh} -t ${ngpus} 

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

    output:
    path("*.ar")
    path("*.png")

    script:
    """
    python3 ${params.fold_script} -i ${peasoup_xml_out} -t pulsarx -p ${params.pulsarx_fold_template} -b ${BEAM} -threads ${psrfold_fil_threads} -ncands ${no_cands_to_fold}
    """

}


workflow {

    // Selecting filterbank_files, target, beam_num, utc_start
    filterbank_channel_with_metadata = Channel
    .fromPath("${params.filterbank_list}")
    .splitCsv(header: true, sep:',') // Ensure the separator is specified if not comma
    .map { row ->
        // Now each row is a Groovy map with column names as keys
        def filterbank_files = row.filterbank_files.trim() // Trim leading and trailing spaces
        def target = row.target.trim()
        def beam_num = row.beam_num.trim()
        def utc_start = row.utc_start.trim().replace(" ", "-") // Replace space with underscore in utc_start

        return tuple(filterbank_files, target, beam_num, utc_start)
    }

    if (params.use_filtool_apsuse == 1) {
        processed_filterbank = filtool_apsuse(filterbank_channel_with_metadata, params.filtool_rfi_filter, params.filtool_threads, params.telescope)
    } else if (params.use_digifil_apsuse == 1) {
        processed_filterbank = digifil_apsuse(filterbank_channel_with_metadata, params.nbits, params.digifil_threads, params.digifil_decimate_time_apsuse)
    }

    // Check if neither processing flags are set to 1
    if (params.use_filtool_apsuse != 1 && params.use_digifil_apsuse != 1) {
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

        if (params.APSUSE_SEARCH == 1) {
            nearest_two_output = nearest_power_of_two_calculator_apsuse(updated_filterbank_channel)
            peasoup_output = peasoup_apsuse(nearest_two_output, params.dm_file, params.total_cands_limit, params.min_snr, params.acc_start, params.acc_end, params.ram_limit_gb, params.nh, params.ngpus, params.kill_file)
            pulsarx_output = fold_peasoup_cands_pulsarx(peasoup_output, params.psrfold_fil_threads, params.no_cands_to_fold)

            // phase_predictor_output = create_phase_predictor_apsuse(peasoup_output, params.period_start, params.period_end)

            // if (params.parallel_fold == 1) {
            //     restructured_phase_predictor_output = phase_predictor_output.flatMap { fil_file, target_name, beam_name, utc_start, fft_size, xml_file, phase_predictors ->
            //         // Check if phase_predictors is a list (or tuple)
            //         if (phase_predictors instanceof List || phase_predictors.getClass().isArray()) {
            //             phase_predictors.collect { phase_predictor ->
            //                 return [fil_file, target_name, beam_name, utc_start, fft_size, xml_file, phase_predictor]
            //             }
            //         } else {
            //             // If it's a single element, return it as it is
            //             return [[fil_file, target_name, beam_name, utc_start, fft_size, xml_file, phase_predictors]]
            //         }
            //     }
                
            //     apsuse_folds = apsuse_fold_phase_predictor_parallel(restructured_phase_predictor_output, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)
            // }
            // //Folding candidates in serial. Use this when you have only a few spin period cands.
            // else{
            //     apsuse_folds = apsuse_fold_phase_predictor_serial(phase_predictor_output, params.dspsr_apsuse_threads, params.telescope, params.dspsr_apsuse_subint_length, params.dspsr_apsuse_bins)
            // }
            

            // if (params.use_clfd == 1) {
            //     clfd_output = clfd_apsuse_predictor(apsuse_folds, params.target_name, params.qmask, params.qspike, params.clfd_processes)
            //     pdmp_output = pdmp_apsuse_predictor(clfd_output, params.target_name, params.nchan, params.nsubint, params.nbins)
            // }
            // else {
            // pdmp_output = pdmp_apsuse_predictor(apsuse_folds, params.target_name, params.nchan, params.nsubint, params.nbins)
            // }

       
        }
        
        }
        
    
        
       
    

    
        
        


 
    



     



