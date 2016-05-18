library(seqinr)

# set up directory structure:
projectname = "benchmark"
chromosome = "chr21"

# use following homeDir for R shell:
#homeDir = "/home/jamtor"
# use following homeDir for R Studio:
homeDir = "/Users/jamestorpy/clusterHome"
projectDir = paste0(homeDir, "/projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
inDir = paste0(resultsDir, chromosome)

setwd(inDir)
"The inDir is:"
inDir

inFile = paste0(chromosome, ".fa")
"The inFile is:"
inFile

#chr = read.fasta(inFile, as.string = T)
#chrvec = unlist(chr)

#chr = read.table(inFile, header = F)

splitchr = split(chrvec, ceiling(seq_along(chrvec)/10000))


v = c(1:100)
result = split(v, ceiling(seq_along(v)/10))
result