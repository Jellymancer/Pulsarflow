# mmgps_pulsar_followup
**A Nextflow Pipeline for Follow-Up of Pulsar Discoveries in the MPIfR MeerKAT Galactic Plane Survey (MMGPS)**

## Installation

### Prerequisites
1. **Java and Nextflow Installation**: Ensure Java is installed on your system. Then, install Nextflow by following the instructions at [Nextflow Getting Started Guide](https://www.nextflow.io/docs/latest/getstarted.html).

### Setting Up the Repository
2. **Environment File**: Inside the repository, create a `.env` file. Populate it with the MMGPS Database credentials (username, password, host, etc.).
3. **Configuration**: Edit `nextflow.config`. Specify `target_name` (e.g., J0737-3039), `beam_name` (e.g. cfbf00000), the spin period range in seconds, and the UTC start and end times. If no UTC times are provided, the software will analyze all available observations.
4. Update the DM range to search with peasoup in the file `dm_range_peasoup.sh`.
5. **Singularity Images**: Update `nextflow.config` with the paths to your Singularity images. You can use the default ones if you are a user of the APSUSE cluster.


### Usage Modes
6. **Operating Modes**: This software supports four modes - APSUSE SEARCH/FOLD MODE and/or PTUSE SEARCH/FOLD Mode. They can be used independently or concurrently.
7. **Initial Pulsar Search**: Without an existing Ephemeris for your pulsar, set `APSUSE_SEARCH=1`, or `PTUSE_SEARCH=1` and `APSUSE_EPH_FOLD=0`, `PTUSE_EPH_FOLD=0`. This configuration will run `filtool`, `peasoup`, and automatically fold all candidates within your specified spin period range in the `nextflow.config` file using `dspsr`.
8. **Folding with Preliminary Ephemeris**: After generating TOAs and obtaining a preliminary ephemeris for APSUSE/PTUSE observations, set `APSUSE_EPH_FOLD` and/or `PTUSE_EPH_FOLD` to `1`. This will fold all candidates with your ephemeris. This might also be useful when you iteratively re-fold your observations as your ephemeris improves over time.
9. **Folding Known Pulsars**: To fold known pulsars in a pointing, update `ephemeris_files_dir` with the path to your directory containing ephemeris files for known pulsars. Set `APSUSE_EPH_FOLD=1` or `PTUSE_EPH_FOLD=1`.

## Running the Pipeline

- **Locally**: Run the pipeline locally using the command:

```
nextflow run main.nf -profile local

```


- **HTCondor Cluster**: To run on the HTCondor cluster, use:

```
nextflow run main.nf -profile condor

```

