# Moore Lab 16S DADA2 Tutorial
# November 18, 2024
# Last Updated: ANB on 11/18/24

# Setup your environment----
# This tutorial assumes you have downloaded the appropriate practice data
# You can read more about fastq files here: 
# https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211

## Load libraries----
library(dada2); library(dplyr); library(ShortRead); library(patchwork)

# Run the following if you need to install dada2
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("dada2")

# Read and inspect your reads----
# Set your path
path="/Users/ashleybulseco/Dropbox/microbial-ecology-DBS/bioinformatics-tutorial/moore-group-tutorials/practice-sequences"
list.files(path)

# Need to make sure this matches what your sequences look like
# Add the .gz if your files are still gunzipped
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))

# To make sure the forward and reverse reads are the same
count_reads <- function(file) {
  fq <- readFastq(file)
  length(fq)
} # count reads function

# Requires the package "ShortRead"
counts_before <- data.frame(
  sample = basename(fnFs),
  forward_reads = sapply(fnFs, count_reads),
  reverse_reads = sapply(fnRs, count_reads)
)

write.csv(counts_before, "read_counts_for_16S.csv")

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
list(sample.names)

# Inspect read quality profiles
P1 <- plotQualityProfile(fnFs[1:11]) # change to your number of samples, max
P2 <- plotQualityProfile(fnRs[1:11])
P1
P2

# Filter & Trim----
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Trimming primers
FWD_PRIMER_LEN <-19
REV_PRIMER_LEN <-20

# Here, you are trimming your reads for quality
# Removing primers
# Truncating the ends of the reads
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     trimLeft = c(FWD_PRIMER_LEN, REV_PRIMER_LEN),
                     truncLen = c(150, 150),  # Example values
                     maxN = 0,  # Discard reads with any Ns
                     maxEE = c(2, 2),  # Allow max 2 expected errors
                     truncQ = 2,  # Truncate reads at the first quality score <= 2
                     compress = TRUE,  # Compress the output files
                     verbose = TRUE, 
                     multithread = TRUE) 
head(out)
saveRDS(out,"out.rds")
# out <- readRDS("out.rds")

# Learn the error rates----
# DADA 2's algorithm is based on modeling the error rates in the sequencing process. 
# It leverages a parametric error model (err), which is different for each individual 
# amplicon dataset. The learnErrors method learns this error model from the data, by
# alternating estimation of the error rates and inference of sample composition until 
# they converge on a jointly consistent solution.
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
saveRDS(errF,"errF.rds")
saveRDS(errR,"errR.rds")
# errF <-readRDS("errF.rds")
# errR <-readRDS("errR.rds")

# Visualize the errors
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# Dereplicate----
# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
saveRDS(derepFs, "derepFS.rds")
saveRDS(derepRs, "derepRs.rds")
# To read back in
# derepFs <- readRDS("derepFs.rds")
# derepRs <- readRDS("derepRs.rds")
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Run the DADA2 algorithm----
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
saveRDS(dadaFs, "dadaFs.rds")
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaRs, "dadaRs.rds")
dadaFs[[1]]
# To read back in
# dadaFs <- readRDS("dadaFs.rds")
# dadaRs <- readRDS("dadaRs.rds")

# Merge your forward and reverse reads----
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])
saveRDS(mergers, "mergers.rds")
# To read back in
# mergers <- readRDS("mergers.rds")

# Construct a sequence table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))

# Remove chimeras----
# Chimeras are common but shouldn't be more than 3-4% 
# If high, then remove primers and redo analysis
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab) # This means that ~2.5% were chimeric
saveRDS(seqtab.nochim, "seqtab-nochim.rds")
# seqtab.nochim <- readRDS("seqtab-nochim.rds")

# Do a quick check in----
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assigning taxonomy----
# Be sure to put the database in your working directory or another known location
# You can download proper database here: https://benjjneb.github.io/dada2/training.html
taxa <- assignTaxonomy(seqtab.nochim, "~/tax/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)
saveRDS(taxa, "taxa-table.rds")
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
