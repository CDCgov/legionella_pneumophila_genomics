#!/bin/bash



###############################################################################
#### Helper Functions ####
###############################################################################

## ****************************************************************************
## Usage description should match command line arguments defined below
usage () {
    echo "Usage: $(basename "$0")"
    echo "  --docker => Docker Image"
    echo "  --r1 => Forward Reads"
    echo "  --r2 => Reverse Reads"
    echo "  --isolate => Isolate to Process"
    echo "  --output => Output Directory"
    echo "  --help => Display this help message"
}
## ****************************************************************************



###############################################################################
## SCRIPT_DIR: directory of current script
###############################################################################
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
## ****************************************************************************



###############################################################################
#### Parse Command-Line Arguments ####
###############################################################################

getopt --test > /dev/null
if [ $? -ne 4 ]; then
    echo "`getopt --test` failed in this environment."
    exit 1
fi

## ****************************************************************************
## Command line options should match usage description
OPTIONS=
LONGOPTIONS=help,docker:,r1:,r2:,isolate:,output:,
## ****************************************************************************

# -temporarily store output to be able to check for errors
# -e.g. use "--options" parameter by name to activate quoting/enhanced mode
# -pass arguments only via   -- "$@"   to separate them correctly
PARSED=$(\
    getopt --options=$OPTIONS --longoptions=$LONGOPTIONS --name "$0" -- "$@"\
)
if [ $? -ne 0 ]; then
    # e.g. $? == 1
    #  then getopt has complained about wrong arguments to stdout
    usage
    exit 2
fi

# read getopt's output this way to handle the quoting right:
eval set -- "$PARSED"

## ****************************************************************************
## Set any defaults for command line options
DOCKER=lp_serogroup:0.1
## ****************************************************************************

## ****************************************************************************
## Handle each command line option. Lower-case variables, e.g., ${file}, only
## exist if they are set as environment variables before script execution.
## Environment variables are used by Agave. If the environment variable is not
## set, the Upper-case variable, e.g., ${FILE}, is assigned from the command
## line parameter.
while true; do
    case "$1" in
        --help)
            usage
            exit 0
            ;;
        --docker)
            DOCKER=$2
            shift 2
            ;;
        --r1)
            R1=$2
            shift 2
            ;;
        --r2)
            R2=$2
            shift 2
            ;;
        --isolate)
            ISOLATE=$2
            shift 2
            ;;
        --output)
            OUTPUT=$2
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option"
            usage
            exit 3
            ;;
    esac
done
## ****************************************************************************

## ****************************************************************************
## Log any variables passed as inputs
echo "Docker: ${DOCKER}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Isolate: ${ISOLATE}"
echo "Output: ${OUTPUT}"
## ****************************************************************************



###############################################################################
#### Validate and Set Variables ####
###############################################################################

## ****************************************************************************
## Add app-specific logic for handling and parsing inputs and parameters

# DOCKER parameter

if [ -z "${DOCKER}" ]; then
    echo "Docker required"
    echo
    usage
    exit 1
fi

# R1 input

if [ -z "${R1}" ]; then
    echo "R1 required"
    echo
    usage
    exit 1
fi
R1_FULL=$(readlink -f "${R1}")
R1_DIR=$(dirname "${R1_FULL}")
R1_BASE=$(basename "${R1_FULL}")

# R2 input

if [ -z "${R2}" ]; then
    echo "R2 required"
    echo
    usage
    exit 1
fi
R2_FULL=$(readlink -f "${R2}")
R2_DIR=$(dirname "${R2_FULL}")
R2_BASE=$(basename "${R2_FULL}")

# ISOLATE parameter

if [ -z "${ISOLATE}" ]; then
    echo "Isolate required"
    echo
    usage
    exit 1
fi

# OUTPUT parameter

if [ -z "${OUTPUT}" ]; then
    echo "Output required"
    echo
    usage
    exit 1
fi
OUTPUT_FULL=$(readlink -f "${OUTPUT}")
OUTPUT_DIR=$(dirname "${OUTPUT_FULL}")
OUTPUT_BASE=$(basename "${OUTPUT_FULL}")

## ****************************************************************************



###############################################################################
#### App Execution Preparation ####
###############################################################################

## ****************************************************************************
## Add logic to prepare environment for execution
mkdir -p ${OUTPUT_FULL}
ulimit -s 5248800
## ****************************************************************************



###############################################################################
#### App Execution ####
###############################################################################

## ****************************************************************************
## Add logic to execute app
docker run -v ${R1_DIR}:/data1 -v ${R2_DIR}:/data2 -v ${OUTPUT_FULL}:/output --rm --privileged ${DOCKER} --r1=/data1/${R1_BASE} --r2=/data2/${R2_BASE} --isolate=${ISOLATE} 2>&1 | tee ${OUTPUT_FULL}/pipeline-${ISOLATE}.log
## ****************************************************************************

