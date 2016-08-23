### reads_mapped_vs_read_length ###

# This script takes the percentage of reads mapped uniquely from artificially generated
# libraries of different read length and plots them as % reads mapped vs read length

# load packages:
library(ggplot2)

# set up directory structure:
projectname = "benchmark"

# if fetching all files from project (not just one sample type), leave 'samplename' blank (i.e. "")
samplename = "long.artfastqgen.reads"
inType = "bwa."

# use following homeDir for R shell:
#homeDir = "/home/jamtor"
# use following homeDir for R Studio:
homeDir = "/Users/jamestorpy/clusterHome"
projectDir = paste0(homeDir, "/projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
inDir = paste0(resultsDir, inType, samplename)
inDir
setwd(inDir)

# fetch unique mapper percentage files, unique IDs, numeric IDs for ordering purposes,
# and number of files:
inFiles = grep("REMOVEME", list.files(path=inDir, pattern = "*unique_mappers.txt",recursive = T), inv=T, value=T)
#inFiles = grep("4900", inFiles, inv=T, value=T)
numericIDs=as.numeric(sub("bp.unique_mappers.txt", "", basename(inFiles)))
no_files=length(inFiles)

# create a dataframe with uniqueIDs of inFiles in column 1, 2nd column empty and 
#prepared for numeric values, column 3 as numeric values of read lengths for sorting:
df=data.frame(numericIDs, column2 = numeric(no_files))

# create vector of uniqueIDs for mapper percentage files and put it in column 2 of df:
i=1
for (file in inFiles) {
  print(file)
  percentage=readChar(file, (file.info(file)$size - 4))
  print(percentage)
  df[i,2]=percentage
  i=i+1
}

# order df by numericIDs and remove numericIDs column afterwards:
df=df[order(df$numericIDs),]

# name columns
colnames(df) = c("read_length", "uniquely_mapped_reads")

df$read_length = as.factor(df$read_length)
df$uniquely_mapped_reads = as.numeric(df$uniquely_mapped_reads)

pdf(paste0(inDir, "/artfastqgen.bwa.length_vs_unique_reads.pdf"))
plot = ggplot(data=df, aes(x=read_length, y=uniquely_mapped_reads, group = 1)
) + geom_point(
) + xlab("Read length (bp)"
) + coord_cartesian(ylim=c(50, 100)
) + ylab("Uniquely mapped reads (%)"
) + scale_x_discrete(breaks=c(50, 75, 120, 200, 300, 400, 700, 900, 2000, 5000, 7500)    
) + scale_y_continuous(breaks=c(50, 60, 70, 80, 90, 1000), labels = c(50, 60, 70, 80, 90, 1000)
) + geom_line()
plot
dev.off()

save.image("mapped_vs_length_computed.RData")

