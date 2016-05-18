#!/bin/bash

module load gi/bwa/0.7.12

numcores=6

#make directory hierachy
projectname="benchmark"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome directory
genome="hg38_ercc"

genomeDir="$homeDir/genomes/$genome"
genomeFile="$genomeDir/$genome.fa"

#input/output types
samplenames=( "artfastqgen.reads" )
outType="bwa"

#extension of files to be used:
inExt=".fastq"

#scripts/logs directory
scriptsPath="$projectDir/scripts/lr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo "The logDir is:"
echo $logDir

#get in/outPaths for all samples:
for samplename in ${samplenames[@]}; do

	#input/output:
	inPath="$resultsDir/$samplename"
	outPath="$resultsDir/$outType.$samplename"
	
	echo -e
	echo This is the inPath:
	echo $inPath

#fetch file names of all projects and put into an array:
	i=0
	files=( $(ls $inPath/2000bp30x$inExt) )
	for file in ${files[@]}; do
		echo -e
		echo The file used is: $file
		filesTotal[i]=$file
		let i++;
	done;

#fetch the inFiles and create an outDir based on their uniqueID:
	j=0
	echo -e
	echo Total files = ${#filesTotal[@]}
	while [ $j -lt ${#filesTotal[@]} ]; do
		inFile=${filesTotal[$j]}
		uniqueID=`basename $inFile | sed s/$inExt//`
		outDir=$outPath/$uniqueID
			
		mkdir -p $outDir

		echo -e
		echo This is the uniqueID:
		echo $uniqueID
		echo -e
		echo This is the outDir:
		echo $outDir

#align reads of input file with BWA, output into .bam files:
      	bwa_line="bwa mem -t $numcores $genomeFile $inFile > $outDir/$uniqueID.sam"

      	echo -e
      	echo "This is the bwa_line:"
        echo $bwa_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
        qsub -N bwa$uniqueID -hold_jid TRIMGALORE_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $bwa_line
			j=$(($j+1))

	done;
done;
