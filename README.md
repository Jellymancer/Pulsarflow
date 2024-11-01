# Pulsarflow
**A Nextflow Pipeline for Pulsar Searches using Peasoup**

## Installation

### Prerequisites
1. **Java and Nextflow Installation**: Ensure Java is installed on your system. Then, install Nextflow by following the instructions at [Nextflow Getting Started Guide](https://www.nextflow.io/docs/latest/getstarted.html).

### Setting Up the Repository
2. **Configuration**: Edit `nextflow.config`. The values for all of the steps can be set therein. The filterbank files to be searched and some auxilleraly values are present in the targets.txt file.
3. **Singularity Images**: Update `nextflow.config` with the paths to your Singularity images. All required images can be pulled from `https://hub.docker.com/u/jellymancer`. 


### Usage Modes
4. **Segmented search**: The script will split the given filterbanks into N chunks before searching.
5. **DM splits**: By adjusting the DM splits values, the search can be run for a number of separate DM chunks. This is usefull when the DM range is large as peasoup tends to have memory issues when searching large DM values.
6. **Search**: The filterbanks (split or whole) are searched with peasoup.
7. **Sifing**: Each pointing is sifted, removing possible RFI and minimizing to no. of times that an identical candidate is folded in different beams.
8. **Candidate Folding**: The script will automatically fold N top SNR candidates (specified in `nextflow.config`) for each beam.
9. **PICS scoring and candyjar preparation**: Finally, the script will score all of the folded candidates and create a candyjar ready tar archive. Check out `https://github.com/vivekvenkris/CandyJar`.

## Running the Pipeline

- **Locally**: Run the pipeline locally using the command:

```
nextflow run main.nf -profile local

```


- **HTCondor Cluster**: To run on the HTCondor cluster, use:

```
WARNING. Currently not supported!
nextflow run main.nf -profile condor

```

- **SLURM Cluster**: To run on the HTCondor cluster, use:

```
nextflow run main.nf -profile slurm

```

