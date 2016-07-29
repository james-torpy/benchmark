#!/bin/bash

numcores=1

#make directory hierachy
projectname="benchmark"
samplename="reads.art_illumina"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output directories
inDir="$resultsDir/$samplename"
outType="trimgalore"
outPath="$resultsDir/$samplename.$outType"

#scripts/log directories
scriptsPath="$projectDir/scripts/sr_QC_map_count"
logDir="$scriptsPath/logs"
mkdir -p $logDir

#extension of files to be used
inExt=".fq.gz"


#fetch file names

i=0

#make 'files' = all files in the input folder with the extension specified as an array
#list all files in the array
#make filesTotal' = number of files processed
#let 1 be added to i each time a file goes through this process as a counter of how many files have been
#processed and a way of assigning each a number:
files=( $(ls $inDir/*$inExt) )
for file in ${files[@]} ;do
	echo -e
	echo The file used is: $file
	echo -e
	filesTotal[i]=$file;
	let i++;
done;

#print how many files there are:
j=0
echo -e
echo The total number of files is: ${#filesTotal[@]}
echo -e


#set up conditions to perform trimgalore analysis on paired raw files
#fetch directory names for files in pairs
#the following function is carried through while the number of files processed 'j' is less than the total
#number of files:

while [ $j -lt ${#filesTotal[@]} ]; do

	inFile1=${files[$j]}

#set up directories specific to each file pair being analysed:
### Need to account for when second pair file appears in folder before first - this sed line doesn't work in this case! ###
	uniqueID=`basename $inFile1 | sed s/_1$inExt//g`
	outDir=$outPath/$uniqueID/


	inFile2=`echo $inFile1 | sed s/_1$inExt/_2$inExt/g`

	mkdir -p $outDir
	
	echo -e
	echo This is the uniqueID:
	echo $uniqueID
	echo -e
	echo This is inFile1:
	echo $inFile1
	echo -e
	echo This is inFile2:
	echo $inFile2
	echo -e
	echo This is the logDir:
	echo $logDir
	echo -e
	echo This is the outDir: $outDir

#echo the command to process the files with trim galore/fastqc:
	trimgalore_line="trim_galore $inFile1 $inFile2 --fastqc --paired --retain_unpaired --length 16 -o $outDir"
	echo -e
	echo The trimgalore_line is:
	echo $trimgalore_line
	echo -e

#submit the job to the cluster as a binary file with name trimgalore_$samplename, creating a log in $logDir
#which includes reported errors:
	qsub -N TRIMGALORE_$uniqueID -wd $logDir -b y -j y -R y -P GenomeInformatics -pe smp $numcores -V $trimgalore_line

	j=$(($j+2))

done;

