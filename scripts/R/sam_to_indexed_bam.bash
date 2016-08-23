#!/bin/bash

### sam_to_bam ###

# This script takes sam files from BWA output of artificial reads and
# converts them to bam

numcores=1

# specify the number of lines chromosome had been split into:
split_no="16390_lines"

#make directory hierachy
projectname="benchmark"

homeDir="/home/jamtor"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output types
samplename="long.artfastqgen.reads"
	inType="bwa"

#extension of files to be used:
inExt=".sam"
outExt=".bam"

#input/output:
inDir="$resultsDir/$inType.$samplename/from_$split_no"_chunks
	
echo -e
echo "This is the inDir:"
echo $inDir

#scripts/logs directory
scriptsPath="$projectDir/scripts/R"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "This is the logDir:"
echo $logDir

#	fetch file names of all projects, put them into an array,
# and convert them into bam files:
files=( $(ls $inDir/**/*bp$inExt) )
for inFile in ${files[@]}; do

	outFile=`echo $inFile | sed "s/$inExt/$outExt/g"`
	outPrefix=`echo $inFile | sed "s/$inExt//g"`
	uniqueID=`basename $inFile | sed "s/$inExt//g"`

	echo -e
	echo "The inFile is:"
	echo $inFile
	echo -e
	echo "The outFile is:"
	echo $outFile

	viewline="samtools view -b -S $inFile > $outFile"
	sortline="samtools sort -T $outPrefix.sorted -o $outPrefix.sorted.bam $outFile"
	indexline="index $outPrefix.sorted.bam"

	echo -e
	echo "The viewline is:"
	echo $viewline
	echo -e
	echo "The sortline is:"
	echo $sortline
	echo -e
	echo "The indexline is:"
	echo $indexline
	echo -e

	qsub -N v_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $viewline
	qsub -N s_$uniqueID -hold_jid v_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $sortline
	qsub -N i_$uniqueID -hold_jid s_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $indexline
	
done;

uniqueFiles=( $(ls $inDir/**/*unique$inExt) )
for inFile in ${uniqueFiles[@]}; do

	uniqueID=`basename $inFile | sed "s/$inExt//g"`
	uniquePrefix=`echo $inFile | sed "s/$inExt//g"`
	originalFile=`echo $inFile | sed "s/.unique//g"`
	outFile=`echo $inFile | sed "s/$inExt/$outExt/g"`
	outPrefix=`echo $inFile | sed "s/$inExt//g"`
	

	echo -e
	echo "The inFile is:"
	echo $inFile
	echo -e
	echo $uniqueID
	echo $uniquePrefix
	echo -e
	echo "The original samfile is:"
	echo $originalFile
	echo -e
	echo "The outFile is:"
	echo $outFile

	headerline="cat $originalFile | head -3 >> $uniquePrefix"\_wheader.sam
	catline="cat $inFile >> $uniquePrefix"\_wheader.sam
	viewline2="samtools view -b -S $uniquePrefix"_wheader".sam > $outFile"
	sortline2="samtools sort -T $uniquePrefix.sorted -o $outPrefix.sorted.bam $outFile"
	indexline2="samtools index $outPrefix.sorted.bam"
	
	echo -e
	echo "The headerline is:"
	echo $headerline
	echo -e
	echo "The catline is:"
	echo $catline
	echo -e
	echo "The viewline is:"
	echo $viewline2
	echo -e
	echo "The sortline is:"
	echo $sortline2
	echo -e
	echo "The indexline is:"
	echo $indexline2
	echo -e

	qsub -N h_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $headerline
	qsub -N c_$uniqueID -hold_jid h_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $catline
	qsub -N v_$uniqueID -hold_jid c_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $viewline2
	qsub -N s_$uniqueID -hold_jid v_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $sortline2
	qsub -N i_$uniqueID -hold_jid s_$uniqueID -b y -wd $logDir -j y -R y -P GenomeInformatics -pe smp $numcores -V $indexline2
	


done;