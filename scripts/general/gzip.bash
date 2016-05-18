#!/bin/bash

#make directory hierachy
projectname="benchmark"
samplename="reads"
inExt=".art_illumina"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"
inPath="$resultsDir/$samplename$inExt"

#create an array of the input .fq files:
inFiles=( $(ls $inPath/*.fq) )
echo -e
echo These are the inFiles:
echo ${inFiles[@]}

#scripts/logs directory
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"
mkdir -p $logDir

#fetch the input .fq files, make a copy of each to retain a gunzipped file appended by .tmp,
#and gzip the original input .fq files:
for file in ${inFiles[@]}; do
	echo -e
	echo The file used is: $file

	copy_line="cp $file "$file.tmp""

	echo -e
	echo This is the copy line:
	echo $copy_line

	gzip_line="gzip $file"

	echo -e
	echo This is the gzip line:
	echo $gzip_line

#rename the .tmp files to the original .fq file names:
	rename_line="mv "$file.tmp" $file"

	echo -e
	echo This is the rename_line:
	echo $rename_line

#define uniqueIDs for each sample to name the cluster jobs:
    uniqueID=`basename $file`

#submit each job to the cluster, holding each job until the previous corresponding jobs are
#complete:
	#qsub -N COPY_$uniqueID -b y -wd $logDir -j y -R y -pe smp 1 -V $copy_line

	qsub -N GZIP_$uniqueID -hold_jid COPY_$uniqueID -b y -wd $logDir -j y -R y -pe smp 1 -V $gzip_line

	#qsub -N RENAME_$uniqueID -hold_jid GZIP_$uniqueID -b y -wd $logDir -j y -R y -pe smp 1 -V $rename_line

done;
