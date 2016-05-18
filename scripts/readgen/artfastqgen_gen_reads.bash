#!/bin/bash

numcores=6

projectname="benchmark"
genomename="hg38_ercc"

#specify the chromosome to generate reads from, the peak coverage mean over each region,
#and the read lengths to generate for each output:
chromosome="chr21"
coverage=5
read_lengths=( 2000 4900 )

#make directory hierachy:
homeDir="/home/jamtor"
binDir="$homeDir/local/bin"
genomeDir="$homeDir/genomes/$genomename"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output directories:
artfastqgen="$binDir/ArtificialFastqGenerator.jar"
genomeFile="$genomeDir/$genomename.fa"
outDir="$resultsDir/reads.artfastqgen"

mkdir -p $outDir

#scripts/log directories
scriptsPath="$projectDir/scripts/readgen"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo -e
echo "The artfastqgen file is:"
echo $artfastqgen
echo -e
echo "The genomeFile is:"
echo $genomeFile
echo -e
echo "The outDir is:"
echo $outDir
echo -e
echo "The logDir is:"
echo $logDir

#output chromosome sequence into .txt file
#index the genomeFile using samtools:
indexline="samtools faidx $genomeFile"

#extract chromosome sequence from genomeFile:
#extractline="samtools faidx $genomeFile $chromosome>$outDir/$chromosome.fa"
extractline="samtools faidx $genomeFile $chromosome>$outDir/$chromosome.fa"

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

#until [ -f $outDir/$chromosome.txt ]
#do
#     sleep 180
#done

#echo -e
#echo "Chromosome file found"
###
###

#count the length in bp of the chromosomes included in the range:
#chr_length=`wc -m $outDir/$chromosome.txt`

#echo -e
#echo "The chr_length is:"
#echo $chr_length

inFile="$resultsDir/$chromosome.fa"

echo -e
echo "This is the inFile:"
echo $inFile

for length in ${read_lengths[@]} ;do
	
	echo -e
	echo "The read_length used is:"
	echo $length

	outPrefix="$outDir/"$length"bp"

	readgenline="java -jar $artfastqgen -O $outPrefix -R $inFile -S '>' -F1 true -F2 true -CMP $coverage -CMPGC 0.4 -CSD 0.2 -GCC false -GCR 10000 -N 10000 \
	-RCNF 2 -RL $length -TLM 9000 -TLSD 60"
	
	echo -e
	echo "The readgenline is:"
	echo $readgenline
	echo -e

	qsub -N readgen$length -hold_jid extractchr -wd $logDir -P GenomeInformatics -b y -j y -R y -pe smp $numcores -V $readgenline
done;
