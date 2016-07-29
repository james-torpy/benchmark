#!/bin/bash

#This script generates error-free simulated reads from a reference fasta file. Read lengths
#are specified by the user in a vector. A specific chromosome is chosen by the user to
#generate the reference fasta file.

numcores=30

projectname="benchmark"
genomename="hg38_ercc"

#specify the chromosome to generate reads from, the peak coverage mean over each region,
#the read lengths to generate for each output:
chromosome="chr21"
coverage=30
read_lengths=( 50 75 120 200 )

#make directory hierachy:
homeDir="/home/jamtor"
binDir="$homeDir/local/bin"
genomeDir="$homeDir/genomes/$genomename"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output directories:
genomeFile="$genomeDir/$genomename.fa"
outDir="$resultsDir/err_free.art_illumina.reads"

mkdir -p $outDir

#scripts/log directories
scriptsPath="$projectDir/scripts/readgen"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "The genomeFile is:"
echo $genomeFile
echo -e
echo "The outDir is:"
echo $outDir
echo -e
echo "The logDir is:"
echo $logDir

#output chromosome sequence into .fa file
#index the genomeFile using samtools:
indexline="samtools faidx $genomeFile"

#extract chromosome sequence from genomeFile:
extractline="samtools faidx $genomeFile $chromosome>$resultsDir/$chromosome.fa"

echo -e
echo "The indexline is:"
echo $indexline
echo -e
echo "The extractline is:"
echo $extractline
echo -e

#submit above lines to the cluster:
#qsub -N indexgenome -wd $logDir -b y -j y -R y -pe smp 1 -V $indexline
#qsub -N extractchr -hold_jid indexgenome -wd $logDir -b y -j y -R y -pe smp 1 -V $extractline

###work out better way of doing the following###
###also make script skip this if chromosome file already exists###
#wait for chromosome file to be created:
#echo "Holding until chromosome file is created"

#until [ -f $outDir/$chromosome.fa ]
#do
#     sleep 180
#done

#echo -e
#echo "Chromosome file found"
###
###

inFile="$resultsDir/$chromosome.fa"

echo -e
echo "This is the inFile:"
echo $inFile

for length in ${read_lengths[@]} ;do
	
	echo -e
	echo "The read_length used is:"
	echo $length

	outPrefix="$outDir/"$length"bp"

	readgenline="art_illumina -i $inFile -l $length -f $coverage -p -m 250 -s 20 -ef -o $outPrefix"

	echo -e
	echo "The readgenline is:"
	echo $readgenline
	echo -e

	#qsub -N readgen$length -hold_jid extractchr -wd $logDir -P GenomeInformatics -b y -j y -R y -pe smp $numcores -V $readgenline
done;
