#!/bin/bash

### artfastqgen_gen_reads_two_thousand_and_above.bash ###

# This script takes a chromosome split into fasta files of different lengths and generates 
# reads of specified lengths from these. It is split for processing purposes and
# can be rejoined using the following scripts: rename_artificial_reads.bash, 
# join_artificial_reads.bash. For reads 900bp or less artfastqgen_gen_reads.bash is
# recommended as processing is not an issue andchromosomes do not need to be split.

numcores=6

projectname="benchmark"
genomename="hg38_ercc"

# specify the number of lines chromosome has been split into:
split_no="16390_lines"

# specify the chromosome to generate reads from, the peak coverage mean over each region,
# and the read lengths to generate for each output:
chromosome="chr21"
coverage=30
read_lengths=( 50 75 120 200 300 400 700 900 2000 5000 7500 )

#make directory hierachy:
homeDir="/home/jamtor"
binDir="$homeDir/local/bin"
genomeDir="$homeDir/genomes/$genomename"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"
inDir="$resultsDir/$chromosome/$split_no"

#input/output directories:
artfastqgen="$binDir/ArtificialFastqGenerator.jar"
outPath="$resultsDir/long.artfastqgen.reads/from_$split_no"_chunks

mkdir -p $outPath

#scripts/log directories
scriptsPath="$projectDir/scripts/readgen"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "The artfastqgen file is:"
echo $artfastqgen
echo -e
echo "The outDir is:"
echo $outDir
echo -e
echo "The logDir is:"
echo $logDir

# make function to generate reads with each read length specified
readgen () {
	for length in "$@"; do

		echo -e
		echo "The read_length used is:"
		echo $length

		outDir="$outPath/$length"bp

		mkdir -p $outDir

		outPrefix="$outDir/"$uniqueID$length"bp"

		echo -e
		echo "The outPrefix is:"
		echo $outPrefix

		# S = <startSequenceIdentifier> The sequence identifier in the reference after which read generation should begin (must be specified)
		# F1/2 = <fastq1ForQualityScores> First/second fastq file to use for real quality scores, (must be specified if useRealQualityScores = true)
		### CMP = <coverageMeanPeak> The peak coverage mean for a region (default = 37.7).
		# GCC = <GCcontentBasedCoverage> Whether nucleobase coverage is biased by GC content (default = true).
		# GCR = <GCcontentRegionSize> Region size in nucleobases for which to calculate GC content, (default = 150)
		# N = <nucleobaseBufferSize> The number of reference sequence nucleobases to buffer in memory, (default = 5000)
		# RCNF = <readsContainingNfilter> Filter out no "N-containing" reads (0), "all-N" reads (1), "at-least-1-N" reads (2), (default = 0)
		# RL = <readLength> The length of each read, (default = 76).
		# TLM, <templateLengthMean> The mean DNA template length, (default = 210)

		readgenline="java -jar $artfastqgen -O $outPrefix -R $inFile -S \
		'>' -F1 true -F2 true -CMP $coverage -GCC false -GCR 10000 -N 10000 \
		-RCNF 2 -RL $length -TLM 9000"
	
		echo -e
		echo "The readgenline is:"
		echo $readgenline
		echo -e

		# submit job to cluster:
		qsub -N $uniqueID$length -hold_jid extractchr -wd $logDir -P GenomeInformatics -b y -j y -R y -pe smp $numcores -V $readgenline
	done;
}

# fetch input files and apply above function:
files=`ls $inDir/$chromosome*.fa | grep -v "$inDir/$chromosome.fa"`
for inFile in $files; do

	# create uniqueID for each file to label cluster job:
	uniqueID=`basename $inFile | sed 's/.fa//g'`

	echo -e
	echo "This is the uniqueID:"
	echo $uniqueID

	readgen ${read_lengths[@]}
done;