process filtool {
    label 'filtool'
    container "${params.pulsarx_singularity_image}"

    publishDir "out/CLEAN/${filterbank_channel_with_metadata[2].trim()}/", pattern: "*.fil", mode: 'symlink'

    input:
    val(filterbank_channel_with_metadata)
    val(rfi_filter_list)
    val(rfi_channel_mask_pairs)
    val threads
    val telescope

    output:
    tuple val(filterbank_channel_with_metadata), path("*.fil")
    
    script:
    def outputFile = "${filterbank_channel_with_metadata[1].trim()}_${filterbank_channel_with_metadata[4].trim()}_${filterbank_channel_with_metadata[2].trim()}_${filterbank_channel_with_metadata[3].trim()}_clean"
    def inputFile = "${filterbank_channel_with_metadata[0].trim()}"
    def source_name = "${filterbank_channel_with_metadata[1].trim()}"

    // Only create zap_string if rfi_channel_mask_pairs is not empty
    def zap_string = rfi_channel_mask_pairs ? rfi_channel_mask_pairs.collect { pair -> "--rfi zap ${pair[0]} ${pair[1]}" }.join(' ') : ""

    // Only create rfi_string if rfi_filter_list is not empty
    def rfi_string = rfi_filter_list ? rfi_filter_list.collect { filter -> "--rfi ${filter}" }.join(' ') : ""

    """
    filtool -t ${threads} --telescope ${telescope} ${zap_string} ${rfi_string} -o ${outputFile} -f ${inputFile} -s ${source_name}
    """
}



process nearest_power_of_two_calculator {
    label 'nearest_power_two'
    container "${params.presto5_singularity_image}"

    input:
    tuple path(fil_file), val(split_id), val(target_name), val(pointing_id), val(beam_name), val(utc_start)

    output:
    tuple path(fil_file), val(split_id), val(target_name), val(pointing_id), val(beam_name), val(utc_start), env(nearest_power_of_2)

    script:
    """
    #!/bin/bash
   
    output=\$(readfile ${fil_file})
    echo "\$output"

    value=\$(echo "\$output" | grep "Spectra per file" | awk '{print \$5}')

    log2=\$(echo "l(\$value)/l(2)" | bc -l)
    decimal_part=\$(echo "\$log2" | awk -F"." '{print "0."\$2}')
    rounded_log2=\$(echo "\$log2" | awk -F"." '{if (\$2 >= 35) print \$1+1; else print \$1}')

    nearest_power_of_2=\$((2**\$rounded_log2))
    echo \$nearest_power_of_2 > nearest_power_of_two.txt

    
    """
}


process calculate_max_acceleration{
    label 'max_acceleration'
    container "${params.presto5_singularity_image}"

    input:
    tuple path(fil_file), val(split_id), val(target_name), val(pointing_id), val(beam_name), val(utc_start), val(nearest_power_of_2)

    output:
    tuple path(fil_file), val(split_id), val(target_name), val(pointing_id), val(beam_name), val(utc_start), val(nearest_power_of_2) ,env(min_accel), env(max_accel)

    """
    output=\$(readfile ${fil_file})
    echo "\$output"

    tobs=\$(echo "\$output" | grep "Time per file (sec)" | awk '{print \$6}')
    max_pobs=\$(echo "scale=2; \$tobs * 10" | bc)

    max_accel=\$(python3 ${params.max_accel_script} -P \${max_pobs})
    min_accel=-\${max_accel}
    """
    
}


process generateDMFiles {

    label 'generate_dm_files'
    container "${params.utility_singularity_image}"

    // publishDir { "${params.dm_out_dir}" }, pattern: "*.txt", mode: 'copy'

    output:
    file("*.txt")

    script:
    """
    python3 ${params.dm_create_script} -ds ${params.dm_start} -de ${params.dm_end} -n ${params.no_dm_splits} -s ${params.dm_step} -d ${params.dm_decimals}
    """
}

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

    publishDir { "out/SEARCH/${utc}/${target_name}/${pointing_id}/${beam_name}/${dm_file.baseName}" }, pattern: "*.xml", mode: 'copy'
    input:
    tuple path(fil_file), val(split_id), val(target_name), val(pointing_id), val(beam_name), val(utc), val(fft_size), val(acc_start), val(acc_end), path(dm_file), val(dm_range)
    val(total_cands_limit)
    val(min_snr)
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
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_file), path(cand_file), path(meta_file)

    output:
    tuple val(target_name), val(pointing_id), val(split_id), val(dm_range), val(beam_name), val(utc), path(fil_file), path("*.csv"), path(meta_file)

    script:
    """
    python3 ${params.sift_script} -c ${params.sift_default_config}  -m ${meta_file} ${cand_file}
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
    container "${params.python2_singularity_image}"
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
