# Legionella pneumophila in silico Serogroup Prediction

Version: 0.0.2

This project is contains a bioinformatics workflow and docker image to analyze Legionella pneumophila whole genome sequencing data to predict serogroup from short read sequences.

## Requirements

In order to execute this workflow please ensure that your computing environment meents the following requirements:

    1. Linux Operating System

    2. SquashFS, required for executing Singularity containers - Most standard Linux distributions have SquashFS installed and enabled by default. However, in the event that SquashFS is not enabled, we recommend that you check with your system administrator to install it. Alternatively, you can enable it by following these instructions (warning: these docs are for advanced users): https://www.tldp.org/HOWTO/html_single/SquashFS-HOWTO/

    3. Docker

## Running the Pipeline Script

```
ulimit -s 5248800
./pipeline.sh --reference=./supportFiles/Phila_NC_002942.fna --gff=./supportFiles/NC_002942.gff --r1=./test-samples/sample_R1.fastq --r2=./test-samples/sample_R2.fastq --isolate=sample --output=./output
```
## Building the Docker Container

### CDC Users

For users in the CDC environment, build the docker container by running the `build.sh` script:

```
./build.sh --tag=lp_serogroup:0.1
```

The `--tag` parameter is optional and has a default value of `lp_serogroup:0.1`.

### External Users

For external users, build the container using ``Dockerfile.external``:

```
docker build -t lp_serogroup:0.1 -f Dockerfile.external .
```

## Docker UID and GID

The docker container will create the output files with the same UID and GID as the output directory, so it is important to create the output directory as yourself (not as root) prior to running the docker container. The wrapper script automatically creates the output directory. 

## Running the Docker Container

The container can be run directly as follows:

```
mkdir output   # the output directory must be created prior to running the docker container
ulimit -s 5248800
docker run -v $(pwd)/test-samples:/data -v $(pwd)/output:/output --privileged --rm lp_serogroup:0.1 --r1=/data/sample_R1.fastq --r2=/data/sample_R2.fastq --isolate=sample
```

Or the container can be run from the wrapper script:

```
./wrapper.sh --docker=lp_serogroup:0.1 --r1=./test-samples/sample_R1.fastq --r2=./test-samples/sample_R2.fastq --isolate=sample --output=./output
```

The wrapper script creates the output directory and sets the ulimit.


### Developed by

[Shatavia Morrison](https://github.com/SMorrison42)


[John Phan](https://github.com/jhphan)
