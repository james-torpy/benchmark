#!/bin/bash

### sam_to_bam ###

# This script takes sam files from BWA output of artificial reads and
# converts them to bam

numcores=6

#make directory hierachy
projectname="benchmark"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output types
samplename="artfastqgen.reads"
inType="bwa"

#extension of files to be used:
inExt=".sam"
outExt=".bam"

#input/output:
inDir="$resultsDir/$inType.$samplename"
	
echo -e
echo "This is the inDir:"
echo $inDir

#	fetch file names of all projects, put them into an array,
# and convert them into bam files:
files=( $(ls $inDir/**/*bp$inExt) )
for inFile in ${files[@]}; do

	outFile=`echo $inFile | sed "s/$inExt/$outExt/g"`
	outPrefix=`echo $inFile | sed "s/$inExt//g"`
	uniqueID=`basename $inFile | sed "s/$inExt//g"`

	echo -e
	echo "The file used is:"
	echo $inFile
	echo -e
	echo "The outFile is:"
	echo $outFile
	echo -e
	echo "The outPrefix is:"
	echo $outPrefix
	echo -e
	echo "The uniqueID is:" $uniqueID

	samtools view -b -S $inFile > $outFile
	samtools sort -T $outPrefix.sorted -o $outPrefix.sorted.bam $outFile
	samtools index $outPrefix.sorted.bam

done;