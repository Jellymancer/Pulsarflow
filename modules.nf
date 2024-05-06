process filtool {
    label 'filtool'
    container "${params.pulsarx_singularity_image}"

    input:
    val(filterbank_channel_with_metadata)
    val rfi_filter
    val threads
    val telescope
    val(channel_mask)

    output:
    tuple val(filterbank_channel_with_metadata), path("*.fil")
    

    script:
    def outputFile = "${filterbank_channel_with_metadata[1].trim()}_${filterbank_channel_with_metadata[3].trim()}_${filterbank_channel_with_metadata[2].trim()}_clean"
    def inputFile = "${filterbank_channel_with_metadata[0].trim()}"
    def source_name = "${filterbank_channel_with_metadata[1].trim()}"
    """
    # Get the first file from the input_data string
    first_file=\$(echo ${inputFile} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"

    # Prepare the rfi mask options
    mask_option=""
    if [[ -n "${channel_mask}" ]]; then
        IFS=',' read -ra ADDR <<< "${channel_mask}"
        for i in "\${ADDR[@]}"; do
            IFS=':' read -ra PAIR <<< "\$i"
            mask_option+=" --rfi zap \${PAIR[0]} \${PAIR[1]}"
        done
    fi

    if [[ \${file_extension} == "sf" ]]; then
        filtool -psrfits --scloffs -t ${threads} --telescope ${telescope} \${mask_option} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
    else 
        filtool -t ${threads} --telescope ${telescope} \${mask_option} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
    fi

    """
}


process nearest_power_of_two_calculator {
    label 'nearest_power_two'
    container "${params.pulsarx_singularity_image}"

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc_start)

    output:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc_start), env(nearest_power_of_2)

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
