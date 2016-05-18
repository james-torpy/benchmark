#!/bin/bash

#This script generates pacbio simulated reads from a reference fasta file. Read lengths
#are specified by the user in a vector. A specific chromosome is chosen by the user to
#generate the reference fasta file.

numcores=6

projectname="benchmark"
genomename="hg38_ercc"

#specify the chromosome to generate reads from, the peak coverage mean over each region,
#the read lengths to generate for each output:
chromosome="chr21"
coverage=30
read_lengths=( 200 300 400 700 900 2000 4900 )

#make directory hierachy:
homeDir="/home/jamtor"
genomeDir="$homeDir/genomes/$genomename"
projectDir="$homeDir/projects/$projectname"

#input/output directories:
genomeFile="$genomeDir/$genomename.fa"
outPath="$projectDir/results/pbsim.reads"
error_profile="$homeDir/local/lib/pbsim-1.0.3/data/model_qc_clr"

mkdir -p $outPath

#scripts/log directories
scriptsPath="$projectDir/scripts/readgen"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "The genomeFile is:"
echo $genomeFile
echo -e
echo "The logDir is:"
echo $logDir

#output chromosome sequence into .fa file
#index the genomeFile using samtools:
indexline="samtools faidx $genomeFile"

#extract chromosome sequence from genomeFile:
extractline="samtools faidx $genomeFile $chromosome>$outPath/$chromosome.fa"

echo -e
echo "The indexline is:"
echo $indexline
echo -e
echo "The extractline is:"
echo $extractline

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

inFile="$outPath/$chromosome.fa"

echo -e
echo "This is the inFile:"
echo $inFile

for length in ${read_lengths[@]} ;do
	
	echo -e
	echo "The read_length used is:"
	echo $length

	outDir="$outPath/$length"
	mkdir -p $outDir

	echo -e
	echo "The outDir is:"
	echo $outDir

	outPrefix="$outDir/"$length"bp"

	readgenline="pbsim --data-type CLR \
				--depth 20 \
				--model_qc $error_profile \
				$inFile"

	echo -e
	echo "The readgenline is:"
	echo $readgenline
	echo -e

	#qsub -N pbsim$length -hold_jid extractchr -wd $outDir -P GenomeInformatics -b y -j y -R y -pe smp $numcores -V $readgenline
done;
