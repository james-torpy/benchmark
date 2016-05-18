#!/bin/bash

### split_chr_into_ten_kbp_for_artfastqgen ###

# This script takes a chromosome fasta file and splits 
# it into files of roughly 10 kbp. For the generation
# of reads of 2000 bp or longer by artfastqgen.

projectname="benchmark"

# specify the chromosome to split:
chromosome="chr21"

#make directory hierachy:
homeDir="/home/jamtor"
binDir="$homeDir/local/bin"
genomeDir="$homeDir/genomes/$genomename"
projectDir="$homeDir/projects/$projectname"
resultsDir="$projectDir/results"

#input/output directories:
inDir="$resultsDir/$chromosome"
origDir="$inDir/original/"
inFile="$inDir/$chromosome.fa"

mkdir -p $origDir

echo -e
echo "The inFile is:"
echo $inFile

# split the chromosome file into files of 1639 lines (or
# ~ 10 kpb) length:
#split -l 1639 $inFile chr21
split -l 50000 $inFile chr21
# move the original chromosome file into separate directory
mv $inFile $origDir

files="$inDir/chr21*"
for file in $files; do
	# add .fa to the end of each file:
	nfile=$file.fa
	mv $file $nfile
	# add '>' to the start of each file so artfastqgen
	# can read them:
	seq=`cat $nfile`
	echo ">"$seq > $nfile
	# remove spaces in the file sequences caused by the
	# cat command:
	sed 's/ //g' $nfile > $nfile.temp
	# add newline every 61 bp for normal fasta format:
	fold -w 61 $nfile.temp > $nfile
	# remove the temp files:
	echo $nfile.temp
	rm $nfile.temp
done;