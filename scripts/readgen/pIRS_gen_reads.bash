#!/bin/bash

#This script generates minimal-error simulated reads from a reference fasta file. Read lengths
#are specified by the user in a vector. A specific chromosome is chosen by the user to
#generate the reference fasta file.

numcores=25

projectname="benchmark"
genomename="hg38_ercc"

#specify the chromosome to generate reads from, the peak coverage mean over each region,
#the read lengths to generate for each output:
chromosome="chr21"
coverage=30
read_lengths=( 50 75 120 200 300 400 )

#make directory hierachy:
homeDir="/home/jamtor"
binDir="$homeDir/local/bin"
genomeDir="$homeDir/genomes/$genomename"
projectDir="$homeDir/projects/$projectname"

#input/output directories:
pIRS="$binDir/pirs"
genomeFile="$genomeDir/$genomename.fa"
outDir="$projectDir/raw_files/pIRS"

mkdir -p $outDir

#scripts/log directories
scriptsPath="$projectDir/scripts/readgen"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "The pIRS file is:"
echo $pIRS
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
extractline="samtools faidx $genomeFile $chromosome>$outDir/$chromosome.fa"

echo -e
echo "The indexline is:"
echo $indexline
echo -e
echo "The extractline is:"
echo $extractline
echo -e

#submit above lines to the cluster:
qsub -N indexgenome -wd $logDir -b y -j y -R y -pe smp 1 -V $indexline
qsub -N extractchr -hold_jid indexgenome -wd $logDir -b y -j y -R y -pe smp 1 -V $extractline

###work out better way of doing the following###
###also make script skip this if chromosome file already exists###
#wait for chromosome file to be created:
echo "Holding until chromosome file is created"

until [ -f $outDir/$chromosome.fa ]
do
     sleep 180
done

echo -e
echo "Chromosome file found"
###
###

inFile="$outDir/$chromosome.fa"

echo -e
echo "This is the inFile:"
echo $inFile

for length in ${read_lengths[@]} ;do
	
	echo -e
	echo "The read_length used is:"
	echo $length

	outPrefix="$outDir/"$length"bp"

	readgenline="pirs simulate -i $inFile -l $length -x $coverage -e 0 -a 0 -g 0 -q 0 -M 0 \
	-s /home/jamtor/local/lib/pIRS_111/src/pirs/src/Profiles/Base-Calling_Profiles/HG00702.PE91.matrix.gz -o $outPrefix"

	echo -e
	echo "The readgenline is:"
	echo $readgenline
	echo -e

	#qsub -N readgen$length -hold_jid extractchr -wd $logDir -P GenomeInformatics -b y -j y -R y -pe smp $numcores -V $readgenline
done;
