#!/bin/bash

numcores=6

#make directory hierachy
projectname="benchmark"

# specify the number of lines chromosome had been split into:
split_no="16390_lines"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output types
samplename="long.artfastqgen.reads"
inType="bwa"

#extension of files to be used:
inExt=".sam"

#scripts/logs directory
scriptsPath="$projectDir/scripts/lr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "This is the logDir:"
echo $logDir

#input/output:
inPath="$resultsDir/$inType.$samplename/from_$split_no"_chunks
	
echo -e
echo "This is the inPath:"
echo $inPath

#fetch file names of all projects and put into an array:
i=0
files=( $(ls $inPath/**/*$inExt) )
for inFile in ${files[@]}; do

	echo -e
	echo The file used is: $inFile

	uniquePrefix=`echo $inFile | sed "s/$inExt//"`

	echo -e
	echo "This is the uniquePrefix:"
	echo $uniquePrefix

	samtools view -Sq 0 $inFile > $uniquePrefix.clean.sam
	samtools view -Sq 1 $inFile > $uniquePrefix.unique.sam
	total_no=`wc -l $uniquePrefix.clean.sam | awk '{print $1;}'`
	unique_no=`wc -l $uniquePrefix.unique.sam | awk '{print $1;}'`
	unique_mappers=`bc <<< "scale = 4; ($unique_no / $total_no)*100"`
	echo $unique_mappers"%" > $uniquePrefix\_unique_mappers.txt
done;

