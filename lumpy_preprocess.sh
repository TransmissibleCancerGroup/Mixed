#!/bin/bash -e

############################################################
#  Program: lumpyexpress
#  Author: Colby Chiang (cc2qe@virginia.edu)
############################################################
set -eo pipefail

# source the paths to the binaries used in the script
source_binaries() {
    if [[ -e $1 ]]
    then
	echo "Sourcing executables from $1 ..."
	if [[ $1 == /* ]]
	then
	    source $1
	else
	    source ./$1
	fi
    else
	echo "Config file $1 not found. Attempting to auto-source executables"
    # general
    LUMPY_HOME=$(dirname $(readlink -f $(which lumpyexpress)))
    if [ ! -d $LUMPY_HOME/scripts ]; then
        LUMPY_HOME=$LUMPY_HOME/..
    fi
    echo $LUMPY_HOME

	LUMPY=`which lumpy || true`
	SAMBLASTER=`which samblaster || true`
	SAMBAMBA=`which sambamba || true`
	SAMTOOLS=`which samtools || true`
	# python 2.7 or newer, must have pysam, numpy installed
	PYTHON=`which python || true`

        # python scripts
	PAIREND_DISTRO=$LUMPY_HOME/scripts/pairend_distro.py
	BAMGROUPREADS=$LUMPY_HOME/scripts/bamkit/bamgroupreads.py
	BAMFILTERRG=$LUMPY_HOME/scripts/bamkit/bamfilterrg.py
	BAMLIBS=$LUMPY_HOME/scripts/bamkit/bamlibs.py
    fi
}

# ensure that the require python modules are installed before
# beginning analysis
check_python_modules() {
    PYTHON_TEST=$1
    echo -e "\nChecking for required python modules ($PYTHON_TEST)..."

    $PYTHON_TEST -c "import imp; imp.find_module('pysam')"
    $PYTHON_TEST -c "import imp; imp.find_module('numpy')"
}

## usage
usage() {
    echo "
usage:   lumpy_preprocess [options]

options:
     -B FILE  full BAM file(s) (comma separated) (required)
     -h       show this message
"
}

# set defaults
LUMPY_DIR=`dirname $0`
CONFIG="$LUMPY_DIR/lumpyexpress.config"
THREADS=1
ANNOTATE=0
MIN_SAMPLE_WEIGHT=4
TRIM_THRES=0
EXCLUDE_BED=
TEMP_DIR=""
GENOTYPE=0
READDEPTH=0
VERBOSE=0
KEEP=0
OUTPUT=""
MAX_SPLIT_COUNT=2
MIN_NON_OVERLAP=20
PROB_CURVE=""
SPL_BAM_STRING=""
DISC_BAM_STRING=""
DEPTH_BED_STRING=""
VERBOSE=1

while getopts ":hB:" OPTION
do
    case "${OPTION}" in
	h)
	    usage
	    exit 0
	    ;;
	B)
	    FULL_BAM_STRING="$OPTARG"
	    ;;
    esac
done

# parse the BAM strings
FULL_BAM_LIST=($(echo $FULL_BAM_STRING | tr "," " "))
SPL_BAM_LIST=($(echo $SPL_BAM_STRING | tr "," " "))
DISC_BAM_LIST=($(echo $DISC_BAM_STRING | tr "," " "))
DEPTH_BED_LIST=($(echo $DEPTH_BED_STRING | tr "," " "))

OPTIND=0

# Check the for the relevant binaries
source_binaries $CONFIG

if [[ -z "$LUMPY" ]]
then
    usage
    echo -e "Error: lumpy executable not found. Please set path in $LUMPY_DIR/lumpyexpress.config file\n"
    exit 1
elif [[ -z  "$PAIREND_DISTRO" ]]
then
    usage
    echo -e "Error: pairend_distro.py executable not found. Please set path in $LUMPY_DIR/lumpyexpress.config file\n"
    exit 1
elif [[ -z "$BAMFILTERRG" ]]
then
    usage
    echo -e "Error: bamfilterrg.py executable not found. Please set path in $LUMPY_DIR/lumpyexpress.config file\n"
    exit 1
fi

# $SAMT will be either sambamba or samtools, depending on which is available
if [[ ! -z "$SAMBAMBA" ]]
then
    SAMT="$SAMBAMBA"
    SAMT_STREAM="$SAMBAMBA view -f bam -l 0"
    SAMTOBAM="$SAMBAMBA view -S -f bam -l 0"
    SAMSORT="$SAMBAMBA sort -m 1G --tmpdir "
elif [[ ! -z "$SAMTOOLS" ]]
then
    SAMT="$SAMTOOLS"
    SAMT_STREAM="$SAMTOOLS view -u"
    SAMTOBAM="$SAMTOOLS view -S -u"
    SAMSORT="$SAMTOOLS sort -m 1G -T "
else
    usage
    echo -e "Error: neither samtools nor sambamba were found. Please set path of one of these in $LUMPY_DIR/lumpyexpress.config file\n"
    exit 1
fi

# check for required python modules (pysam, numpy)
check_python_modules $PYTHON

# Check that the required files exist
if [[ ${#FULL_BAM_LIST[@]} -eq 0 ]]
then
    usage
    echo -e "Error: -B is required\n"
    exit 1
fi

set +o nounset
for TEST_BAM in ${FULL_BAM_LIST[@]} ${SPL_BAM_LIST[@]} ${DISC_BAM_LIST[@]}
do
    if [[ ! -f $TEST_BAM ]]
    then
	usage
	echo -e "Error: file $TEST_BAM not found.\n"
	exit 1
    fi
done

for TEST_BED in ${DEPTH_BED_LIST[@]}
do
	if [[ -z $(echo "$TEST_BED" | grep ":") ]]
	then
		usage
		echo -e "Error: must specify depths as sample_id:bedpe"
		exit 1;
	fi
	bpath=$(echo "$TEST_BED" | perl -pe 's/^.+://')
	if [[ ! -f "$bpath" ]]; then
		usage
		echo -e "Error: depth bed does not exist: $bpath"
		exit 1
	fi
done
set -o nounset

# default OUTPUT if not provided
if test -z "$OUTPUT"
then
    OUTPUT=`basename "${FULL_BAM_LIST[0]}"`.tmp
fi
OUTBASE=`basename "$OUTPUT"`

# make temporary directory
if [[ $VERBOSE -eq 1 ]]
then
    echo "
    create temporary directory"
fi
if [[ -z $TEMP_DIR ]]
then
    TEMP_DIR=`mktemp -d ${OUTBASE}.XXXXXXXXXXXX`
else
    mkdir -p $TEMP_DIR
fi


cleanup () {
	rm -rf $TEMP_DIR
}
trap cleanup EXIT

####################################
# Extract split and discordant reads
####################################

set +o nounset
# initialize split and discordant bam lists
SPL_BAM_LIST=()
DISC_BAM_LIST=()

# create temp files
mkdir -p $TEMP_DIR/spl $TEMP_DIR/disc

# generate histo files and construct the strings for LUMPY
for i in $( seq 0 $(( ${#FULL_BAM_LIST[@]}-1 )) ); do
    FULL_BAM=${FULL_BAM_LIST[$i]}

    # calc readlength if not provided
    set +o pipefail
    READ_LENGTH=`$SAMT view $FULL_BAM | head -n 10000 | gawk 'BEGIN { MAX_LEN=0 } { LEN=length($10); if (LEN>MAX_LEN) MAX_LEN=LEN } END { print MAX_LEN }'`
    set -o pipefail

    # parse the libraries in the BAM header to extract readgroups from the same library
    LIB_RG_LIST=(`$PYTHON $BAMLIBS $FULL_BAM`)

    # process each library's splitters and discordants
    for j in $( seq 0 $(( ${#LIB_RG_LIST[@]}-1 )) ); do
        SPLITTER=${FULL_BAM%.bam}.spl.sam
        DISCORDS=${FULL_BAM%.bam}.disc.sam

        if [[ "$VERBOSE" -eq 1 ]]; then
            echo -e "$PYTHON $BAMGROUPREADS --fix_flags -i $FULL_BAM -r ${LIB_RG_LIST[$j]} \
    | $SAMBLASTER --acceptDupMarks --excludeDups --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP \
    --splitterFile $SPLITTER --discordantFile $DISCORDS > /dev/null"
            echo -e "$SAMTOBAM $SPLITTER | $SAMSORT $TEMP_DIR/spl -o ${SPLITTER%.sam}.bam /dev/stdin"
            echo -e "$SAMTOBAM $DISCORDS | $SAMSORT $TEMP_DIR/disc -o ${DISCORDS%.sam}.bam /dev/stdin" 
        fi

        $PYTHON $BAMGROUPREADS --fix_flags -i $FULL_BAM -r ${LIB_RG_LIST[$j]} \
        | $SAMBLASTER --acceptDupMarks --excludeDups --addMateTags --maxSplitCount $MAX_SPLIT_COUNT --minNonOverlap $MIN_NON_OVERLAP \
            --splitterFile $SPLITTER --discordantFile $DISCORDS > /dev/null && echo "Samblaster stage success!"

        $SAMTOBAM $SPLITTER | $SAMSORT $TEMP_DIR/spl -o ${SPLITTER%.sam}.bam /dev/stdin && rm $SPLITTER && echo "Split reads sorting stage: success!"
        $SAMTOBAM $DISCORDS | $SAMSORT $TEMP_DIR/disc -o ${DISCORDS%.sam}.bam /dev/stdin && rm $DISCORDS && echo "Discordant reads sorting stage: success!"
        wait
    done
done
exit
