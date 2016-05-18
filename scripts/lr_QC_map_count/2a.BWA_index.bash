#!/bin/bash

module load gi/bwa/0.7.12

#number of cores
numcores=2

#genome directories
genomeName="hg38_ercc"

homeDir="/home/jamtor"
genomeDir="$homeDir/genomes/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"
outDir="$genomeDir/bwa.index"

mkdir -p $outDir

#log directory
projectname="benchmark"

projectDir="$homeDir/projects/$projectname"
scriptsPath="$projectDir/scripts/lr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the outDir
echo $outDir
echo -e

#generate the star reference files:
bwa_index_line="bwa index $genomeFile"

echo This is the bwa_index_line:
echo $bwa_index_line

#submit job with name 'bwaindex' to 15 cluster cores:
qsub -N bwaindex -wd $logDir -b y -cwd -j y -R y -pe smp $numcores -V $bwa_index_line
