#!/bin/bash

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

#scripts/logs directory
scriptsPath="$projectDir/scripts/lr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "This is the logDir:"
echo $logDir

#input/output:
inDir="$resultsDir/$inType.$samplename"
outDir=$inPath
	
echo -e
echo "This is the in/outDir:"
echo $outDir

#fetch file names of all projects and put into an array:
i=0
files=( $(ls $inDir/**/*$inExt) )
for inFile in ${files[@]}; do

	echo -e
	echo The file used is: $inFile

	uniquePath=`echo $inFile | sed 's/$inExt//g'`
	uniqueID=`basename $inFile | sed 's/$inExt//g'`

	echo -e
	echo "This is the uniquePath:"
	echo $uniquePath
	echo -e
	echo "This is the uniqueID:"
	echo $uniqueID

	#samtools view -Sq 0 $inFile > $uniquePath.clean.sam
	#samtools view -Sq 1 $inFile > $uniquePath.unique.sam
	#total_no=`wc -l $uniquePath.clean.sam | sed 's/$uniqueID.clean.sam//g'`
	#unique_no=`wc -l $uniquePath.unique.sam | sed 's/$uniqueID.unique.sam//g'`
	#unique_mappers=`bc <<< "scale = 4; ($unique_no / $total_no)*100"`
	#echo $unique_mappers"%" > $uniqueID_unique_mappers.txt
done;

