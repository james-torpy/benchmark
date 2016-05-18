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

#input/output:
inDir="$resultsDir/$inType.$samplename"
	
echo -e
echo "This is the inDir:"
echo $inDir


#fetch file names of all projects and put into an array:
files=( $(ls $inDir/**/2000bp30x$inExt) )
for inFile in ${files[@]}; do

	echo -e
	echo The file used is: $inFile

	uniquePrefix=`echo $inFile | sed "s/$inExt//g"`
	uniqueID=`basename $inFile | sed "s/$inExt//g"`

	echo -e
	echo "This is the uniquePrefix:"
	echo $uniquePrefix
	echo -e
	echo "This is the uniqueID:"
	echo $uniqueID

	samtools view -Sq 0 $inFile > $uniquePrefix.clean.sam
	samtools view -Sq 1 $inFile > $uniquePrefix.unique.sam
	total_no=`wc -l $uniquePrefix.clean.sam`
	total_no=`echo $total_no | awk '{print $1}'`
	unique_no=`wc -l $uniquePrefix.unique.sam`
	unique_no=`echo $unique_no | awk '{print $1}'`
	unique_mappers=`bc <<< "scale = 4; ($unique_no / $total_no)*100"`
	echo $unique_mappers"%" > $uniquePrefix.unique_mappers.txt
done;

