### extract_dna-type_from_bam.R ###

# This script takes a bam file and extracts specific types of DNA, placing it in a
# separate bam file (e.g. protein coding, lnc, repeat, pseudogene)

# load packages:
library(ggplot2)
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)

# set up directory structure:
projectname <- "benchmark"
annotation <- "hg38_ercc"

# if fetching all files from project (not just one sample type), leave 'samplename' blank (i.e. "")
samplename <- "artfastqgen.reads"
inType <- "bwa."
inExt <- ".sorted.bam"

# use following homeDir for R shell:
#homeDir = "/home/jamtor"
# use following homeDir for R Studio:
homeDir <- "/Users/jamestorpy/clusterHome"
annotDir <- paste0(homeDir, "/genomes/", annotation)
projectDir <- paste0(homeDir, "/projects/", projectname)
resultsDir <- paste0(projectDir, "/results/")
inDir <- paste0(resultsDir, inType, samplename)
inDir
setwd(inDir)

annotFile <- paste0(annotDir, "/gencode_v24_hg38_annotation.gtf")
annotFile

# specify ScanBam parameters:
what <- c("qname", "rname", "strand", "pos", "qwidth", "qual")
flag <- scanBamFlag(isUnmappedQuery=FALSE)
param <- ScanBamParam(what=what, flag=flag)
  
# fetch the bam files:
bamFiles <- grep(list.files(inDir, pattern = paste0("bp", inExt), recursive = T
), pattern = "REMOVEME", inv = T, value = T)
bamFiles <- grep(bamFiles, pattern = ".bai", inv = T, value = T)
print("The bamFiles used are:")
bamFiles

# create annotation in a GRangesList format:
#annot_list <- import(annotFile)
categories <- c("protein_coding", "lincRNA", "pseudogene")
###find out how to get the repeats!###
protein_coding <- reduce(annot_list[annot_list$gene_type=="protein_coding"])
lincRNA <- reduce(annot_list[annot_list$gene_type=="lincRNA"])
pseudogene <- reduce(annot_list[annot_list$gene_type=="pseudogene"])

annotation <- GRangesList()
annotation[["protein_coding"]] <- protein_coding
annotation[["lincRNA"]] <- lincRNA
annotation[["pseudogene"]] <- pseudogene

# load in bam files to ScanBam - the count for each of the reads is in the name of the
# read, after "c":
sampleName <- gsub(".genome.sorted.bam", "", bamFiles)

for(bamFile in bamFiles) {
  paste0("Converting ", bamFile, " to GRanges object")
  split_filename <- str_split(bamFile, "/")
  uniqueID <- split_filename[[1]][1]
  print(paste0("The uniqueID is: ", uniqueID))
  system.time(bam <- scanBam(bamFile, param=param)[[1]])
  
  # convert bam file to GRanges object
  gr <- GRanges(
    seqnames = bam$rname,
    ranges = IRanges(start=bam$pos, width=bam$qwidth),
    strand = bam$strand,
    id = bam$qname
  )
  
  # find overlaps between annotations and bam files:
  sampleCounts <- list()
  for(type in names(annotation)) {
    mat <- findOverlaps(gr, annotation[[type]])
    shortGR <- gr[unique(queryHits(mat))]
    uniqIDs <- shortGR$id
    sampleCounts[[type]] <- length(uniqIDs)
    gr_temp <- gr[!(gr$id %in% uniqIDs)]
    assign(paste0(uniqueID, "_", type, "_gr"), gr_temp)
  }
}

#Check only mapped to one chromosome:
#length(grep(gr, pattern = "chr13", value = T))
#length(grep(gr, pattern = "chr14", value = T))
#length(grep(gr, pattern = "chr15", value = T))
#length(grep(gr, pattern = "chr16", value = T))
#length(grep(gr, pattern = "chr17", value = T))
#length(grep(gr, pattern = "chr18", value = T))
#length(grep(gr, pattern = "chr19", value = T))
#length(grep(gr, pattern = "chr20", value = T))
#length(grep(gr, pattern = "chr21", value = T))
#search <- grep(gr, pattern = "chr21", inv = T, value = T)
#length(grep(search, pattern = "chr13", inv = T, value = T))
