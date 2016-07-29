#!/bin/bash

module load gi/samtools/1.2

numcores=6

#make directory hierachy
projectname="benchmark"
samplename=""

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#genome directory
genome="hg38_ercc"

genomeDir="$homeDir/genomes/$genome"

echo -e
echo This is the genomeDir:
echo $genomeDir
echo -e

#input/output types
samplenames=( "art_illumina.reads" )
inType="trimgalore"
outType="star"

#extension of files to be used:
inExt=".fq.gz"

#scripts/logs directory
scriptsPath="$projectDir/scripts/sr_QC_map_count"
logDir="$scriptsPath/logs"

mkdir -p $logDir

#get in/outPaths for all samples:
for samplename in ${samplenames[@]}; do

	#input/output:
	inPath="$resultsDir/$samplename"
	outPath="$resultsDir/$outType.$samplename"
	
	echo This is the inPath:
	echo $inPath
	echo -e

#fetch file names of all projects and put into an array:
	i=0
	files=( $(ls $inPath/*$inExt | grep -v unpaired) )
	for file in ${files[@]}; do
		echo The file used is: $file
		echo -e
		filesTotal[i]=$file
		let i++;
	done;

#fetch the inFiles and create an outDir based on their uniqueID:
	j=0
	echo Total files = ${#filesTotal[@]}
	echo -e
	while [ $j -lt ${#filesTotal[@]} ]; do
		inFile1=${filesTotal[$j]}
		inFile2=${filesTotal[$(($j+1))]}
		uniqueID=`basename $inFile1 | sed s/_1$inExt//`
		outDir=$outPath/$uniqueID/
			
		mkdir -p $outDir

		echo -e
		echo This is the uniqueID:
		echo $uniqueID
		echo -e
		echo This is the outDir:
		echo $outDir
		echo -e

#align reads of input files with STAR, output into .bam files:
		starJobName="star."$uniqueID
        bamJobName="bam."$uniqueID
        sortJobName="sort."$uniqueID

        indexJobName="index."$uniqueID
        indexStatsJobName="indexstats."$uniqueID
        outSam=$outDir"Aligned.out.sam"
        outBam=$outDir"$uniqueID.bam"
        outSortedBam=$outDir"$uniqueID.sorted.bam"

      	star_line="STAR --runMode alignReads \
       	--readFilesCommand zcat \
     	--genomeDir $genomeDir \
    	--outFilterType BySJout \
      	--outSAMattributes NH HI AS NM MD\
      	--outFilterMultimapNmax 20 \
      	--outFilterMismatchNmax 999 \
      	--outFilterMismatchNoverReadLmax 0.04 \
      	--alignIntronMin 20 \
      	--alignIntronMax 1500000 \
      	--alignMatesGapMax 1500000 \
      	--alignSJoverhangMin 6 \
      	--alignSJDBoverhangMin 1 \
      	--readFilesIn $inFile1 $inFile2 \
      	--outFileNamePrefix $outDir \
      	--runThreadN $numcores \
	--quantMode TranscriptomeSAM \
      	--outFilterMatchNmin 76"

        echo $star_line

#submit jobs to the cluster, creating a log in $logDir which includes reported errors:                	
      qsub -N STAR_$uniqueID -hold_jid TRIMGALORE_$UniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $star_line

			j=$(($j+2))

	done;
done;
