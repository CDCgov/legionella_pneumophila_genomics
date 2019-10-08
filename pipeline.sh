#!/bin/bash



###############################################################################
#### Helper Functions ####
###############################################################################

## ****************************************************************************
## Usage description should match command line arguments defined below
usage () {
    echo "Usage: $(basename "$0")"
    echo "  --reference => Reference Sequence"
    echo "  --gff => Gene Features"
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
LONGOPTIONS=help,reference:,gff:,r1:,r2:,isolate:,output:,
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
#
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
        --reference)
            REFERENCE=$2
            shift 2
            ;;
        --gff)
            GFF=$2
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
echo "Reference: ${REFERENCE}"
echo "GFF: ${GFF}"
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

# REFERENCE input

if [ -z "${REFERENCE}" ]; then
    echo "Reference required"
    echo
    usage
    exit 1
fi
REFERENCE_FULL=$(readlink -f "${REFERENCE}")
REFERENCE_DIR=$(dirname "${REFERENCE_FULL}")
REFERENCE_BASE=$(basename "${REFERENCE_FULL}")
if [ ! -f "${REFERENCE_FULL}" ]; then
    echo "Reference file not found"
    echo
    usage
    exit 1
fi

# GFF

if [ -z "${GFF}" ]; then
    echo "GFF required"
    echo
    usage
    exit 1
fi
GFF_FULL=$(readlink -f "${GFF}")
GFF_DIR=$(dirname "${GFF_FULL}")
GFF_BASE=$(basename "${GFF_FULL}")
if [ ! -f "${GFF_FULL}" ]; then
    echo "GFF file not found"
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
if [ ! -f "${R1_FULL}" ]; then
    echo "R1 file not found"
    echo
    usage
    exit 1
fi

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
if [ ! -f "${R2_FULL}" ]; then
    echo "R2 file not found"
    echo
    usage
    exit 1
fi

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
mkdir -p ${OUTPUT_FULL}/logs
## ****************************************************************************



###############################################################################
#### App Execution ####
###############################################################################

## ****************************************************************************
## Add logic to execute app

### Make Bowtie Index
mkdir -p ${OUTPUT_FULL}/bowtie2-index

echo '01) bowtie2 index'
singularity -s exec -B ${REFERENCE_DIR}:/data -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/bowtie2:2.2.8--py35_2 bowtie2-build /data/${REFERENCE_BASE} /output/bowtie2-index/${REFERENCE_BASE} > ${OUTPUT_FULL}/logs/01-bowtie2-index.stdout 2> ${OUTPUT_FULL}/logs/01-bowtie2-index.stderr || { echo 'bowtie2 index failed'; exit 1; }

### Map Reads to Reference Legionella pneumophila Philadelphia genome
mkdir -p ${OUTPUT_FULL}/map-sam

echo '02) bowtie2 align'
singularity -s exec -B ${R1_DIR}:/data1 -B ${R2_DIR}:/data2 -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/bowtie2:2.2.8--py35_2 bowtie2 -x /output/bowtie2-index/${REFERENCE_BASE} --very-sensitive-local --no-unal -a -1 /data1/${R1_BASE} -2 /data2/${R2_BASE} -S /output/map-sam/${ISOLATE}.sam > ${OUTPUT_FULL}/logs/02-bowtie2-align.stdout 2> ${OUTPUT_FULL}/logs/02-bowtie2-align.stdout || { echo 'bowtie2 align failed'; exit 1; }

### Convert Sam to Bam File
mkdir -p ${OUTPUT_FULL}/map-bam

echo '03) samtools view'
singularity -s exec -B ${REFERENCE_DIR}:/data -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/samtools:1.3.1--3 samtools view -bS -T /data/${REFERENCE_BASE} /output/map-sam/${ISOLATE}.sam -o /output/map-bam/${ISOLATE}.bam > ${OUTPUT_FULL}/logs/03-samtools-view.stdout 2> ${OUTPUT_FULL}/logs/03-samtools-view.stderr || { echo 'samtools view failed'; exit 1; }

echo '04) samtools sort'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/samtools:1.3.1--3 samtools sort -o /output/map-bam/${ISOLATE}-sorted.bam -O bam /output/map-bam/${ISOLATE}.bam > ${OUTPUT_FULL}/logs/04-samtools-sort.stdout 2> ${OUTPUT_FULL}/logs/04-samtools-sort.stderr || { echo 'samtools sort failed'; exit 1; }

echo '05) samtools index'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/samtools:1.3.1--3 samtools index /output/map-bam/${ISOLATE}-sorted.bam > ${OUTPUT_FULL}/logs/05-samtools-index.stdout 2> ${OUTPUT_FULL}/logs/05-samtools-index.stderr || { echo 'samtools index failed'; exit 1; }

mkdir -p ${OUTPUT_FULL}/variants

### Call SNPs with FreeBayes
echo '06) freebayes'
singularity -s exec -B ${REFERENCE_DIR}:/data -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/freebayes:1.2.0--py27h56106d0_4 freebayes -q 20 -p 1 --min-coverage 20 -F 0.75 -j -f /data/${REFERENCE_BASE} /output/map-bam/${ISOLATE}-sorted.bam -v /output/variants/${ISOLATE}-freebayes-SNP-differences.vcf > ${OUTPUT_FULL}/logs/06-freebayes.stdout 2> ${OUTPUT_FULL}/logs/06-freebayes.stderr || { echo 'freebayes failed'; exit 1; }

### Use VCFTools to create a VCF File
echo '07) vcftools'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/vcftools:0.1.14--5 vcftools --vcf /output/variants/${ISOLATE}-freebayes-SNP-differences.vcf --remove-indels --recode --recode-INFO-all --out /output/variants/${ISOLATE}-freebayes-SNP-differences-NO-INDELS > ${OUTPUT_FULL}/logs/07-vcftools.stdout 2> ${OUTPUT_FULL}/logs/07-vcftools.stderr || { echo 'vcftools failed'; exit 1; }

### Use vcfFilter to remove low quality SNPs
echo '08) vcffilter'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/vcflib:1.0.0_rc2--h56106d0_1 vcffilter -f "QUAL > 1" /output/variants/${ISOLATE}-freebayes-SNP-differences-NO-INDELS.recode.vcf > ${OUTPUT_FULL}/variants/${ISOLATE}-vcffilter.vcf 2> ${OUTPUT_FULL}/logs/08-vcffilter.stderr || { echo 'vcffilter failed'; exit 1; }

### GZip VCF file
echo '09) bgzip'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/htslib:1.9--h47928c2_5 bgzip /output/variants/${ISOLATE}-vcffilter.vcf > ${OUTPUT_FULL}/logs/09-bgzip.stdout 2> ${OUTPUT_FULL}/logs/09-bgzip.stderr || { echo 'bgzip failed'; exit 1; }

### Use Tabix command on Gzipped VCF file
echo '10) tabix'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/htslib:1.9--h47928c2_5 tabix -f -p vcf /output/variants/${ISOLATE}-vcffilter.vcf.gz > ${OUTPUT_FULL}/logs/10-tabix.stdout 2> ${OUTPUT_FULL}/logs/10-tabix.stderr || { echo 'tabix failed'; exit 1; }

mkdir -p ${OUTPUT_FULL}/consensus

### Use VCF Consensus to infer SNPs detected by FreeBayes into Isolate Reference Consensus Sequence
echo '11) vcf-consensus'
cat ${REFERENCE_FULL} | singularity -s exec -B ${OUTPUT_FULL}:/output docker://biocontainers/vcftools:v0.1.15_cv2 vcf-consensus /output/variants/${ISOLATE}-vcffilter.vcf.gz > ${OUTPUT_FULL}/consensus/${ISOLATE}.fasta 2> ${OUTPUT_FULL}/logs/11-vcf-consensus.stderr || { echo 'vcf-consensus failed'; exit 1; }

### Convert 60 line character fasta file to one line fasta file
echo '12) awk one-line fasta'
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' ${OUTPUT_FULL}/consensus/${ISOLATE}.fasta > ${OUTPUT_FULL}/consensus/${ISOLATE}-one-line.fasta 2> ${OUTPUT_FULL}/logs/12-awk-one-line-fasta.stderr || { echo 'awk one-line fasta failed'; exit 1; }

### Use Bedtools to identify the depth coverage at each site in consensus sequence
echo '13) bedtools genomecov'
singularity -s exec -B ${REFERENCE_DIR}:/data -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/bedtools:2.23.0--hdbcaa40_3 bedtools genomecov -bga -split -ibam /output/map-bam/${ISOLATE}-sorted.bam -g /data/${REFERENCE_BASE} > ${OUTPUT_FULL}/consensus/${ISOLATE}-per-site-depth-coverage.txt 2> ${OUTPUT_FULL}/logs/13-bedtools-genomecov.stderr || { echo 'bedtools genomecov failed'; exit 1; }

### Use python script to mask sites below 25x depth coverage
echo '14) mask sites'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://smorrison42/lpserogroup_python:latest python /scripts/maskSites.py /output/consensus/${ISOLATE}-one-line.fasta /output/consensus/${ISOLATE}-per-site-depth-coverage.txt /output/consensus/${ISOLATE}-masked-seq.fasta > ${OUTPUT_FULL}/logs/14-mask-sites.stdout 2> ${OUTPUT_FULL}/logs/14-mask-sites.stderr || { echo 'mask sites failed'; exit 1; }

### Use Bedtools to identify gene sequences and their start/end coordinates
echo '15) bedtools getfasta'
singularity -s exec -B ${GFF_DIR}:/data -B ${OUTPUT_FULL}:/output docker://quay.io/biocontainers/bedtools:2.26.0gx--he860b03_3 bedtools getfasta -fi /output/consensus/${ISOLATE}-masked-seq.fasta -bed /data/${GFF_BASE} -name -fo > ${OUTPUT_FULL}/consensus/${ISOLATE}-genes-LPS-regions.fasta 2> ${OUTPUT_FULL}/logs/15-bedtools-getfasta.stderr || { echo 'bedtools getfasta failed'; exit 1; }

mkdir -p ${OUTPUT_FULL}/features

### Use python script to calculate the coverage percentage for genes within the LPS Biosynthesis region based on L. pneumophila Philadelphia gene coords: below 80% denoted as absent;and make input matrix of SNPs and sites
echo '16) gene presence absence matrix'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://smorrison42/lpserogroup_python:latest python /scripts/genePresenceAbsenceMatrix.py /output/consensus/${ISOLATE}-genes-LPS-regions.fasta ${ISOLATE} /output/features/${ISOLATE}-output-genes-LPS.txt > ${OUTPUT_FULL}/logs/15-gene-presence-absence-matrix.stdout 2> ${OUTPUT_FULL}/logs/15-gene-presence-absence-matrix.stderr || { echo 'gene presence absence matrix failed'; exit 1; }

### Use python script to identify nucleotide sites for RF prediciton: sites were identified from training set removing overall invariant sites, need the RF_sites_20180622.txt file
echo '17) extract sites for prediction'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://smorrison42/lpserogroup_python:latest python /scripts/extractSitesforPrediction.py /output/features/${ISOLATE}-output-genes-LPS.txt /output/features/${ISOLATE}-sites-extract.csv > ${OUTPUT_FULL}/logs/17-extract-sites-for-prediction.stdout 2> ${OUTPUT_FULL}/logs/17-extract-sites-for-prediction.stderr || { echo 'extract sites for prediction failed'; exit 1; }

mkdir -p ${OUTPUT_FULL}/predict

### need to set ulimit to 5248800 to run the model in R
echo '18) r model predict'
ulimit -s 5248800
singularity -s exec -B ${OUTPUT_FULL}:/output docker://smorrison42/lpserogroup_rscripts:latest Rscript --max-ppsize=500000 /scripts/R_modelPredict_version_0.1_20180904Model.R /output/features/${ISOLATE}-sites-extract.csv ${ISOLATE} /output/predict/ > ${OUTPUT_FULL}/logs/18-r-model-predict.stdout 2> ${OUTPUT_FULL}/logs/18-r-model-predict.stderr || { echo 'r model predict failed'; exit 1; }

### Merge prediction result with L.pneumophila training model per serogroup sensitivity and specificity percentage to assist with acceptance or rejecting Lp serogroup prediction
echo '19) merge predict results'
singularity -s exec -B ${OUTPUT_FULL}:/output docker://smorrison42/lpserogroup_python:latest python /scripts/mergePredictResults.py /output/predict/${ISOLATE}-predictResults.txt ${ISOLATE} /output/predict/${ISOLATE}-final-prediction-results.txt > ${OUTPUT_FULL}/logs/19-merge-predict-results.stdout 2> ${OUTPUT_FULL}/logs/19-merge-predict-results.stderr || { echo 'merge predict results failed'; exit 1; }

echo 'Pipeline Complete!'
