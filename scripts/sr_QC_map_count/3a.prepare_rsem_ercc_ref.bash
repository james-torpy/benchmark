#!/bin/bash

numcores=20

module load gi/boost/1.53.0

#directory hierachy
projectname="Grant"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"

#genome/annotation file directories
genomeName="hg38_ercc"
annotationName="gencode.v24.annotation"

genomeDir="$homeDir/genomes/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"
annotationFile="$genomeDir/$annotationName/annotationName.gtf"
outDir="$genomeDir/rsem.ref"

mkdir -p $outDir

outFile="$outDir/rsem_hg38_ercc_ref"

#scripts/log directory
scriptsPath="$projectDir/scripts/sr_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir:
echo $outDir
echo -e

#Set up conditions
rsem_ref_line="rsem-prepare-reference --gtf $annotationFile $genomeFile $outDir"

echo This is the rsem_ref_line:
echo $rsem_ref_line

qsub -N RSEM_ref_$genome -wd $logDir -b y -j y -R y -pe smp $numcores -V $rsem_ref_line
