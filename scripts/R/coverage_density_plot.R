### coverage_density_plot ###

# This script takes a .bam file and generates a density plot of the coverage
# of the regions the reads were mapped to 

# load packages:
library("Rsamtools")
library("BSgenome.Hsapiens.UCSC.hg38")

# set up directory structure:
projectname = "benchmark"

# if fetching all files from project (not just one sample type), leave 'samplename' blank (i.e. "")
samplename = "long.artfastqgen.reads"
inType = "bwa."

# specify how many lines of the reference the chunks reads were generated from:
chunk_lines <- 16390

# use following homeDir for R shell:
#homeDir = "/home/jamtor"
# use following homeDir for R Studio:
homeDir = "/Users/jamestorpy/clusterHome"
projectDir = paste0(homeDir, "/projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
inDir = paste0(resultsDir, inType, samplename, "/from_", chunk_lines, "_lines_chunks")
paste0("The inDir is: ", inDir)
setwd(inDir)

# determine inFiles:
inFiles = list.files(inDir, pattern = ".sorted.bam", recursive = T)
inFiles = grep("bai", inFiles, inv=T, value=T)
inFiles = grep("REMOVEME", inFiles, inv=T, value=T)
inFiles = grep("unique", inFiles, inv=T, value=T)
paste0("The inFiles are: ")
inFiles

# specify column names of bam dataframe
what = c("qname","flag","rname","strand","pos","qwidth")
#define parameters of bam scan:
param = ScanBamParam(what=what)

# create empty vector for bam objects:
bamObjects <-  vector(mode = "character")

i = 1
for (file in inFiles) {
  # check file:
  print(paste0("The file used is: ", file))
  
  # index the bam file:
  indexBam(file)
  
  # create uniqueID to identify each bam file:
  uniqueID <- basename(sub(".sorted.bam", "", file))
  print(paste0("The uniqueID is: ", uniqueID))
  
  # scan the bam file into a list (of 6 attribute lists of 1 list each) with uniqueID:
  bam=scanBam(file,param=param)
  assign(paste0("_bam", uniqueID), bam)
  
  # sanity check:
  paste0("The dimensions of _bam", uniqueID, " are: ")
  dim(paste0("_bam", uniqueID))
  
  # add bam dataframe to bamObjects vector (of lists from each bam file of 6 lists of 1 list each):
  bamObjects[i] <- paste0("_bam", uniqueID)
  i = i + 1
}

### understand the following:

# function for collapsing the list of list of 6 lists into a list of 6 lists:
unlist_it <- function (x){
  x1 <- x[[1L]]
  if(is.factor(x1)){
    structure(unlist(x), class = "factor", levels = levels(x1))
  } else {
    do.call(c, x)
  }
}

###label each thing printed for following loop:

for (object in bamObjects) {
  
  # create a unique ID to identify each bam file:
  uID <- basename(sub("_bam", "", object))
  print(paste0("The uniqueID is: ", uID))
  
  # grab object from vector:
  object <- get(object)
  
  # store names of bam fields:
  bam_field <- names(object[[1]])
  
  # go through each bam field and unlist:
  list <- lapply(bam_field, function(y) unlist_it(lapply(object, "[[", y)))
  
  # store as data frame:
  bam_df <- do.call("DataFrame", list)
  names(bam_df) <- bam_field
  
  # check dimensions of bam_df:
  print(dim(bam_df))
  
  # find out how many reads mapped to negative strand, and assign entries to
  # variable:
  print(table(bam_df$strand == '-'))
  neg_reads <- subset(bam_df, strand == '-')
  neg_reads <- neg_reads$pos
  
  # sanity check:
  print(length(neg_reads))
  
  # find out how many reads mapped to positive strand, and assign the mapped
  # positions to a variable:
  print(table(bam_df$strand == '+'))
  pos_reads <- subset(bam_df, strand == '+')
  pos_reads <- pos_reads$pos
  
  # sanity check:
  print(length(pos_reads))
  
  # additional sanity check:
  NA_reads <- bam_df[is.na(bam_df$strand),]
  nrow(NA_reads)
  print(nrow(bam_df) == length(neg_reads) + length(pos_reads) + nrow(NA_reads))
  
  # assign unique IDs to neg_reads and pos_reads objects:
  assign(paste0(uID, "_neg_reads"), neg_reads)
  assign(paste0(uID, "_pos_reads"), pos_reads)
}

for (object in bamObjects) {
  # create IDs to identify each bam file and read lengths of these:
  uID <- basename(sub("bam_", "", object))
  print(paste0("The uniqueID is: ", uID))
  
  # fetch neg and pos read locations
  negReads <- get(paste0("bam_", uID, "_neg_reads"))
  posReads <- get(paste0("bam_", uID, "_pos_reads"))
  
  # calculate the densities:
  neg_density <- density(negReads)
  pos_density <- density(posReads)

  # change the negative density to a negative scale for plot:
  neg_density$y <- neg_density$y * -1
  
  # generate density plot for both neg and pos mapped reads:
  pdf(paste0(inDir, "/coverage_density_", uID, ".pdf"))
  plot(pos_density,
       ylim = range(c(neg_density$y, pos_density$y)),
       main = paste0("Coverage of mapped ", uID, " artificial reads"),
       xlab = "Chromosome 21",
       col = 'blue',
       lwd = 2.5)
  lines(neg_density, lwd = 2.5, col = 'red')
  dev.off()
}

# save final RData object:
save.image(paste0(projectDir, "/Robjects/coverage_density_plot_final.RData"))
