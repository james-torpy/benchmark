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
