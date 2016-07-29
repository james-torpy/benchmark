samtools view -Sq 0 900bpsample.1.bam > 900bpsample.1.clean.sam
samtools view -Sq 1 900bpsample.1.bam > 900bpsample.1.unique.sam

total_no=`wc -l 900bpsample.1.clean.sam | sed 's/900bpsample.1.clean.sam//g'`
unique_no=`wc -l 900bpsample.1.unique.sam | sed 's/900bpsample.1.unique.sam//g'`

unique_proportion=`echo $(( unique_no/total_no ))`

# > unique_mappers.txt