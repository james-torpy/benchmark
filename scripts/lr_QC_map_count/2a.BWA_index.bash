#!/bin/bash

module load gi/bwa/0.7.12

#number of cores
numcores=2

# chromosome directory:
projectname="benchmark"
chromosome="chr21"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"
chrDir="$resultsDir/$chromosome/original"
chrFile="$chrDir/$chromosome.fa"

# log directory
projectname="benchmark"

projectDir="$homeDir/projects/$projectname"
scriptsPath="$projectDir/scripts/lr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo This is the chrDir:
echo $chrDir
echo -e
echo This is the chrFile:
echo $chrFile
echo -e

#generate the star reference files:
bwa_index_line="bwa index $chrFile"

echo This is the bwa_index_line:
echo $bwa_index_line

#submit job with name 'bwaindex' to 15 cluster cores:
qsub -N bwaindex -wd $logDir -b y -cwd -j y -R y -pe smp $numcores -V $bwa_index_line
