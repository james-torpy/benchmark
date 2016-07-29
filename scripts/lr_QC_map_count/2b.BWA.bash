#!/bin/bash

module load gi/bwa/0.7.12

numcores=6

# specify the number of lines chromosome had been split into to generate the reads:
split_no="16390_lines"

#make directory hierachy
projectname="benchmark"
chromosome="chr21"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"
chrDir="$resultsDir/$chromosome/original"
chrFile="$chrDir/$chromosome.fa"

#input/output types
samplenames=( "long.artfastqgen.reads" )
outType="bwa"

#extension of files to be used:
inExt=".combined.fastq"

#scripts/logs directory
scriptsPath="$projectDir/scripts/lr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "This is the chrFile:"
echo $chrFile
echo -e
echo "This is the logDir:"
echo $logDir

#get in/outPaths for all samples:
for samplename in ${samplenames[@]}; do

	#input/output:
	inPath="$resultsDir/$samplename/from_$split_no"_chunks
	outPath="$resultsDir/$outType.$samplename/from_$split_no"_chunks
	
	echo -e
	echo This is the inPath:
	echo $inPath

#fetch file names of all projects and put into an array:
	i=0
	files=( $(ls $inPath/**/combined/*$inExt) )
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
		inFile1=${filesTotal[$j]}
		inFile2=${filesTotal[$j+1]}
		uniqueID=`basename $inFile1 | sed s/.1$inExt//`
		outDir="$outPath/$uniqueID"\bp
			
		mkdir -p $outDir

		echo -e
		echo This is the uniqueID:
		echo $uniqueID
		echo -e
		echo This is the outDir:
		echo $outDir

#align reads of input file with BWA, output into .bam files:
      	bwa_line="bwa mem -t $numcores $chrFile $inFile1 $inFile2 > $outDir/$uniqueID"bp.sam

      	echo -e
      	echo "This is the bwa_line:"
        echo $bwa_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
        qsub -N bwa$uniqueID -hold_jid TRIMGALORE_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $bwa_line
			j=$(($j+2))

	done;
done;
