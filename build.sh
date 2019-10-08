#!/bin/bash



###############################################################################
#### Helper Functions ####
###############################################################################

## ****************************************************************************
## Usage description should match command line arguments defined below
usage () {
    echo "Usage: $(basename "$0")"
    echo "  --tag => Docker image tag"
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
LONGOPTIONS=help,tag:,
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
TAG=lp_serogroup:0.1
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
        --tag)
            TAG=$2
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
echo "Tag: ${TAG}"
## ****************************************************************************



###############################################################################
#### Validate and Set Variables ####
###############################################################################

## ****************************************************************************
## Add app-specific logic for handling and parsing inputs and parameters

# TAG parameter

if [ -z "${TAG}" ]; then
    echo "Docker image tag required"
    echo
    usage
    exit 1
fi

## ****************************************************************************



###############################################################################
#### App Execution Preparation ####
###############################################################################

## ****************************************************************************
# check if CDC specific certs exist and copy them
rm -rf ./certs
mkdir ./certs
if [ -d "/usr/share/pki/ca-trust-source/anchors/" ]; then
    cp /usr/share/pki/ca-trust-source/anchors/* ./certs
fi
# if no certs, make sure there is at least one file for the docker copy
touch ./certs/dummy.txt
## ****************************************************************************



###############################################################################
#### App Execution ####
###############################################################################

## ****************************************************************************
# build docker image
docker build -t ${TAG} -f Dockerfile.cdc .
## ****************************************************************************


