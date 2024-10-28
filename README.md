# mmgps_pulsar_followup
**A Nextflow Pipeline for Pulsar Searches using Peasoup**

## Installation

### Prerequisites
1. **Java and Nextflow Installation**: Ensure Java is installed on your system. Then, install Nextflow by following the instructions at [Nextflow Getting Started Guide](https://www.nextflow.io/docs/latest/getstarted.html).

### Setting Up the Repository
2. **Configuration**: Edit `nextflow.config`. Specify `target_name` (e.g., J0737-3039), `beam_name` (e.g. cfbf00000), `filterbank_file` (absolute path), the spin period range in seconds, and the UTC start and end times (only used for tagging cleaned filterbank files). Also secify the DM range and the DM subdivisons (if required). T
3. **Singularity Images**: Update `nextflow.config` with the paths to your Singularity images.


### Usage Modes
4. **Initial Pulsar Search**: The script will run `filtool` and `peasoup`. The latter is run  independently for the specified DM subdivions. This can help when memory/cluster time are limited. 
5. **Candidate Folding**: The script will then automatically fold N top SNR candidates (specified in `nextflow.config`) file using `psrfold_fil`.

## Running the Pipeline

- **Locally**: Run the pipeline locally using the command:

```
nextflow run main.nf -profile local

```


- **HTCondor Cluster**: To run on the HTCondor cluster, use:

```
nextflow run main.nf -profile condor

```

- **SLURM Cluster**: To run on the HTCondor cluster, use:

```
nextflow run main.nf -profile slurm

```

