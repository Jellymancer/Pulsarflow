process filtool {
    label 'filtool'
    container "${params.pulsarx_singularity_image}"

    input:
    val(filterbank_channel_with_metadata)
    val rfi_filter
    val threads
    val telescope

    output:
    tuple val(filterbank_channel_with_metadata), path("*.fil")
    

    script:
    def outputFile = "${filterbank_channel_with_metadata[1].trim()}_${filterbank_channel_with_metadata[3].trim()}_${filterbank_channel_with_metadata[2].trim()}"
    def inputFile = "${filterbank_channel_with_metadata[0].trim()}"
    def source_name = "${filterbank_channel_with_metadata[1].trim()}"
    """
    # Get the first file from the input_data string
    first_file=\$(echo ${inputFile} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    if [[ \${file_extension} == "sf" ]]; then
        filtool -psrfits --scloffs -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
    else 
        filtool -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f ${inputFile} -s ${source_name}
    fi

    """
}

process filtool_ptuse {
    label 'filtool'
    container "${params.pulsarx_singularity_image}"

    input:
    tuple val(input_data), val(target_name), val(utc_start)
    val rfi_filter
    val threads
    val telescope

    output:
    tuple val(target_name), val(utc_start), path("*.fil")
    

    script:
    def outputFile = "${target_name}_${utc_start}"
    """
    # Get the first file from the input_data string
    inputFiles=\$(ls -v ${input_data} | tr '\\n' ' ')
    first_file=\$(echo \${inputFiles} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"
    
    if [[ \${file_extension} == "sf" ]]; then
        filtool --psrfits --scloffs -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f \${inputFiles} -s ${target_name}
    else 
        filtool -t ${threads} --telescope ${telescope} -z ${rfi_filter} -o ${outputFile} -f \${inputFiles} -s ${target_name}
    fi

    """
}

process create_phase_predictor {
    label 'create_phase_predictor' 
    container "${params.fold_singularity_image}"

    input:
    tuple path(fil_file), val(target_name), val(beam_name), val(utc_start), val(fft_size), path(xml_file) 
    val(period_start)
    val(period_end)

    output:
    tuple (path(fil_file), val(target_name), val(beam_name), val(utc_start), val(fft_size), path(xml_file), path("predictor_candidate_*"), optional: true)


    script:
    """
    python ${params.phase_predictor} -i ${xml_file} -ps ${period_start} -pe ${period_end} -t ${target_name}
    """
}

process nearest_power_of_two_calculator {
    label 'nearest_power_two'
    container "${params.fold_singularity_image}"

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
process dspsr_fold_phase_predictor {
    label 'fold_phase_predictor' 
    container "${params.fold_singularity_image}"

    input:
    tuple val(fil_file), path(phase_predictor)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)
    val(target_name)

    output:
    path "*.ar"

    script:
    """
    # Extract the basename of the filterbank file
    output_filename=\$(basename ${fil_file} .fil)

    # Extract the basename of the phase_predictor file
    phase_predictor_basename=\$(basename ${phase_predictor})

    # Extract the candidate number from the phase_predictor basename
    candidate_id=\$(echo \${phase_predictor_basename} | sed 's/predictor_candidate_//')

    # Append the candidate_id to the output filename
    output_filename_with_candidate_id=\${output_filename}_candidate_id_\${candidate_id}

    # Run dspsr with the modified output filename
    dspsr -k ${telescope} -t${threads} -P ${phase_predictor} -L ${subint_length} -b${bins} -A -O \${output_filename_with_candidate_id} ${fil_file}
    """
}

process dspsr_fold_phase_predictor_serial {
    label 'fold_phase_predictor' 
    container "${params.fold_singularity_image}"

    input:
    tuple val(fil_file), val(target_name), val(beam_name), val(utc_start), val(fft_size), val(xml_file), val(phase_predictors)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)

    output:
    path "*.ar"

    script:
    """
    #!/bin/bash
    set -ue

    # Correctly preprocess phase_predictors to remove square brackets, spaces, and convert to a Bash array
    predictors_cleaned=\$(echo "${phase_predictors}" | sed 's/\\[//g' | sed 's/\\]//g' | tr -d ' ' | tr ',' '\\n')

    readarray -t phase_predictor_array <<< "\$predictors_cleaned"

    # Function to run dspsr command
    run_dspsr () {
        local fil_file=\$1
        local phase_predictor=\$2
        local telescope=${telescope}
        local threads=${threads}
        local subint_length=${subint_length}
        local bins=${bins}

        output_filename=\$(basename \$fil_file .fil)
        phase_predictor_basename=\$(basename \$phase_predictor)
        candidate_id=\${phase_predictor_basename//predictor_candidate_/}
        output_filename_with_candidate_id=\${output_filename}_candidate_id_\${candidate_id}

        dspsr -k \$telescope -t\$threads -P \$phase_predictor -L \$subint_length -b\$bins -A -O \$output_filename_with_candidate_id \$fil_file
    }

    # Check if phase_predictors is a list or a single file
    if [[ \${#phase_predictor_array[@]} -gt 1 ]]; then
        # It's a list, iterate over each element
        for phase_predictor in "\${phase_predictor_array[@]}"; do
            run_dspsr ${fil_file} \$phase_predictor
        done
    else
        # It's a single file
        run_dspsr ${fil_file} "\${phase_predictor_array[0]}"
    fi
    """
}





// process dspsr_fold_phase_predictor_serial {
//     label 'fold_phase_predictor' 
//     container "${params.fold_singularity_image}"

//     input:
//     tuple val(fil_file), val(target_name), val(beam_name), val(utc_start), val(fft_size), val(xml_file), val(phase_predictors)
//     val(threads)
//     val(telescope)
//     val(subint_length)
//     val(bins)

//     output:
//     path "*.ar"

//     script:
//     phase_predictors.flatten().each { phase_predictor ->
//         """
//         echo ${phase_predictor}
//         output_filename=\$(basename ${fil_file} .fil)
//         phase_predictor_basename=\$(basename ${phase_predictor})
//         candidate_id=\$(echo \${phase_predictor_basename} | sed 's/predictor_candidate_//')
//         output_filename_with_candidate_id=\${output_filename}_candidate_id_\${candidate_id}
//         dspsr -k ${telescope} -t${threads} -P ${phase_predictor} -L ${subint_length} -b${bins} -A -O \${output_filename_with_candidate_id} ${fil_file}
//         """
//     }
// }

process dspsr_fold_phase_predictor_parallel {
    label 'fold_phase_predictor' 
    container "${params.fold_singularity_image}"

    input:
    tuple val(fil_file), val(target_name), val(beam_name), val(utc_start), val(fft_size), val(xml_file), path(phase_predictor)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)

    output:
    path "*.ar"

    script:
    """
    output_filename=\$(basename ${fil_file} .fil)
    phase_predictor_basename=\$(basename ${phase_predictor})
    candidate_id=\$(echo \${phase_predictor_basename} | sed 's/predictor_candidate_//')
    output_filename_with_candidate_id=\${output_filename}_candidate_id_\${candidate_id}
    dspsr -k ${telescope} -t${threads} -P ${phase_predictor} -L ${subint_length} -b${bins} -A -O \${output_filename_with_candidate_id} ${fil_file}
    """
}


process dspsr_fold_ephemeris {
    label 'fold_ephemeris' 
    container "${params.fold_singularity_image}"

    input:
    tuple val(file_path), val(target_name), val(beam_name), val(utc_start), path(ephemeris_file)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)

    output:
    path "*.ar"

    script:
    """
    # Find at least one file that matches the wildcard pattern
    found_file=\$(ls -v ${file_path} 2> /dev/null | head -n 1)

    if [[ -z "\$found_file" ]]; then
        echo "Error: No matching files found."
        exit 1
    fi

    # Extracting basename and filename without extension from ephemeris_file
    ephemeris_basename=\$(basename ${ephemeris_file})
    ephemeris_name_no_ext=\${ephemeris_basename%.*}
    
    # Construct the output file name with ephemeris file suffix
    output_filename="${target_name}_${utc_start}_\${ephemeris_name_no_ext}"

    # Check the file extension of the found file
    file_extension="\${found_file##*.}"
    

    # Run the dspsr command
    if [[ \${file_extension} == "sf" ]]; then
        dspsr -scloffs -k ${telescope} -t${threads} -E ${ephemeris_file} -L ${subint_length} -b${bins} -A -O \${output_filename} \$(ls -v ${file_path})
    else
        dspsr -k ${telescope} -t${threads} -E ${ephemeris_file} -L ${subint_length} -b${bins} -A -O \${output_filename} \$(ls -v ${file_path})
    fi
    """
}


process dspsr_fold_ephemeris_ptuse {
    label 'fold_ephemeris'
    container "${params.fold_singularity_image}"

    input:
    tuple val(file_path), val(target_name), val(beam_name), val(utc_start), path(ephemeris_file)
    val(threads)
    val(telescope)
    val(subint_length)
    val(bins)

    output:
    path "*.ar"

    script:
    """
    # Find at least one file that matches the wildcard pattern
    found_file=\$(ls -v ${file_path} 2> /dev/null | head -n 1)

    if [[ -z "\$found_file" ]]; then
        echo "Error: No matching files found."
        exit 1
    fi

    # Extracting basename and filename without extension from ephemeris_file
    ephemeris_basename=\$(basename ${ephemeris_file})
    ephemeris_name_no_ext=\${ephemeris_basename%.*}

    # Initialize an array to keep track of successfully processed files
    declare -a successful_files

    

    # Loop through all files in the specified path
    for input_file in \$(ls -v ${file_path}); do
        # Extract the filename without its path and extension
        filename=\$(basename "\$input_file")
        filename_no_ext="\${filename%.*}"

        
        output_filename="\${filename_no_ext}.ar"

        # Process each file with dspsr. If dspsr fails for a file, the script continues without exiting or throwing an error.
        if dspsr -scloffs -A -E ${ephemeris_file} -L ${subint_length} -b ${bins} -t ${threads} -k ${telescope} -e ar -O "\${filename_no_ext}" "\$input_file" >& /dev/null; then
            # If dspsr succeeds, add the path of the output .ar file to the list of successful files
            successful_files+=("\$output_filename")
        fi

        
    done

    # Check if there are any successful files to process further
    if [[ \${#successful_files[@]} -gt 0 ]]; then
        # Construct the psradd command using the array of successfully processed files
        psradd -o "${target_name}_${utc_start}_\${ephemeris_name_no_ext}.ar" \${successful_files[@]}
        for file in "\${successful_files[@]}"; do
            rm -rfv "\$file"
        done
    else
        echo "No files were successfully processed by dspsr."
    fi
    """
}




process clfd {
    label 'clfd' 
    container "${params.clfd_singularity_image}"

    input:
    path fold_archive
    val target_name
    val qmask
    val qspike
    val processes

    output:
    path "*.clfd"

    script:
    """
    #!/bin/bash
    set -ue

    # Convert fold_archive to an array, splitting by spaces
    read -a fold_archives_array <<< "${fold_archive}"

    # Check the number of elements in fold_archives_array
    if [[ \${#fold_archives_array[@]} -gt 1 ]]; then
        # Multiple archives, process each one separately
        for archive in "\${fold_archives_array[@]}"; do
            clfd --fmt psrchive \$archive --features std ptp lfamp --qmask ${qmask} --despike --qspike ${qspike} --processes ${processes} -o \${PWD}
        done
    else
        # Single archive
        clfd --fmt psrchive ${fold_archive} --features std ptp lfamp --qmask ${qmask} --despike --qspike ${qspike} --processes ${processes} -o \${PWD}
    fi
    """
}


process pdmp {
    label 'pdmp' 
    container "${params.fold_singularity_image}"

    input:
    path fold_archive
    val target_name
    val nchan
    val nsubint
    val nbin

    output:
    path "*.png"

    script:
    """
    #!/bin/bash
    set -ue

    # Convert fold_archive to an array, splitting by spaces
    read -a fold_archives_array <<< "${fold_archive}"

    # Check the number of elements in fold_archives_array
    if [[ \${#fold_archives_array[@]} -gt 1 ]]; then
        # Multiple archives, process each one separately
        for archive in "\${fold_archives_array[@]}"; do
            output_filename=\$(basename \$archive | sed 's/\\.[^.]*\$//')
            pdmp -mc ${nchan} -ms ${nsubint} -mb ${nbin} -g \${output_filename}.png/PNG \$archive
        done
    else
        # Single archive
        output_filename=\$(basename ${fold_archive} | sed 's/\\.[^.]*\$//')
        pdmp -mc ${nchan} -ms ${nsubint} -mb ${nbin} -g \${output_filename}.png/PNG ${fold_archive}
    fi
    """
}


process digifil {
    label 'digifil'
    container "${params.digifil_singularity_image}"

    input:
    val(filterbank_channel_with_metadata)
    val(nbits)
    val(nthreads)
    val(time_decimate_factor)

    output:
    tuple val(filterbank_channel_with_metadata), path("*.fil")

    script:

    def outputFile = "${filterbank_channel_with_metadata[1].trim()}_${filterbank_channel_with_metadata[3].trim()}_${filterbank_channel_with_metadata[2].trim()}.fil"
    def inputFile = "${filterbank_channel_with_metadata[0].trim()}"

    """
    # Get the first file from the inputFile string
    first_file=\$(echo ${inputFile} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"

    # Determine the -scloffs flag based on file extension
    scloffs_flag=""
    if [[ "\${file_extension}" == "sf" ]]; then
        scloffs_flag="-scloffs"
    fi

    # Conditional execution of digifil based on time_decimate_factor
    if [ ${time_decimate_factor} -eq 1 ]; then
        digifil \${scloffs_flag} -b ${nbits} -threads ${nthreads} -d 1 -o ${outputFile} \$(ls -v ${inputFile})
    else
        digifil \${scloffs_flag} -b ${nbits} -threads ${nthreads} -t ${time_decimate_factor} -d 1 -o ${outputFile} \$(ls -v ${inputFile})
    fi
    """
}

process digifil_ptuse {
    label 'digifil'
    container "${params.digifil_singularity_image}"

    input:
    tuple val(input_data), val(target_name), val(utc_start)
    val(nbits)
    val(nthreads)
    val(time_decimate_factor)

    output:
    tuple val(target_name), val(utc_start), path("*.fil")

    script:

    def outputFile = "${target_name}_${utc_start}"

    """
    # Get the first file from the inputFile string
    inputFiles=\$(ls -v ${input_data} | tr '\\n' ' ')
    first_file=\$(echo \${inputFiles} | awk '{print \$1}')

    # Extract the file extension from the first file
    file_extension="\$(basename "\${first_file}" | sed 's/.*\\.//')"

    # Determine the -scloffs flag based on file extension
    scloffs_flag=""
    if [[ "\${file_extension}" == "sf" ]]; then
        scloffs_flag="-scloffs"
    fi

    # Conditional execution of digifil based on time_decimate_factor
    if [ ${time_decimate_factor} -eq 1 ]; then
        digifil \${scloffs_flag} -b ${nbits} -threads ${nthreads} -d 1 -o ${outputFile} \${inputFiles}
    else
        digifil \${scloffs_flag} -b ${nbits} -threads ${nthreads} -t ${time_decimate_factor} -d 1 -o ${outputFile} \${inputFiles}
    fi
    """
}

