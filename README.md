# Legionella pneumophila in silico Serogroup Prediction

Version: 0.2

This project contains a bioinformatics workflow and docker image to analyze Legionella pneumophila whole genome sequencing data to predict serogroup from short read sequences.

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
./build.sh --tag=lpserogroup_prediction:0.2
```

The `--tag` parameter is optional and has a default value of `lpserogroup_prediction:0.2`.

### External Users

For external users, build the container using ``Dockerfile.external``:

```
docker build -t lpserogroup_prediction:0.2 -f Dockerfile.external .
```

## Docker UID and GID

The docker container will create the output files with the same UID and GID as the output directory, so it is important to create the output directory as yourself (not as root) prior to running the docker container. The wrapper script automatically creates the output directory. 

## Running the Docker Container

The container can be run directly as follows:

```
mkdir output   # the output directory must be created prior to running the docker container
ulimit -s 5248800
docker run -v $(pwd)/test-samples:/data -v $(pwd)/output:/output --privileged --rm lpserogroup_prediction:0.2 --r1=/data/sample_R1.fastq --r2=/data/sample_R2.fastq --isolate=sample
```

Or the container can be run from the wrapper script:

```
./wrapper.sh --docker=lpserogroup_prediction:0.2 --r1=./test-samples/sample_R1.fastq --r2=./test-samples/sample_R2.fastq --isolate=sample --output=./output
```
 
 Or you can pull the container from:
 
 ```
 docker pull smorrison42/lpserogroup_prediction:0.2
```

The wrapper script creates the output directory and sets the ulimit.

## Developed by

[Shatavia Morrison](https://github.com/SMorrison42)


[John Phan](https://github.com/jhphan)

## License

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the Surveillance Platform Disclaimer and Code of Conduct. For more information about CDC's privacy policy, please visit http://www.cdc.gov/privacy.html.

## Contributing

Anyone is encouraged to contribute to the repository by forking and submitting a pull request. (If you are new to GitHub, you might start with a basic tutorial.) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the Apache Software License v2 or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page are subject to the Presidential Records Act and may be archived. Learn more at http://www.cdc.gov/other/privacy.html.

## Records

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the CDC web site.

## Notices

Please refer to CDC's Template Repository for more information about contributing to this repository, public domain notices and disclaimers, and code of conduct.
