
# Make sure your R is the current version

install.packages("installr")
library(installr)
if(!require(installr)) {
  install.packages("installr"); 
  require(installr)
}
updateR()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")

# Intall required packages if you have not already done so. 
# Make sure you have the Bioconductor repository installed

BiocManager::install("dada2")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")

# Update and oad required packages, plus any others you may want to use

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

# Set working path
# Add correct file path that matches your unique path

path <- "C:/Users/af206/Desktop/09082021_AF_DADA2/HT_seqs"  
list.files(path)

#path1 <- "S:/Research/Peter Allen Lab/Ashley Fletcher/Allen Microbiome/Nature_Matters_Arising/Miller_Data/fastq_files/HT"
#list.files(path1)

fnFs <- sort(list.files(path, pattern = "_001.fastq.gz", full.names = TRUE)) #Change file name accordingly
fReads <- file.path(path, "fReads", basename(fnFs))
readFastq(pattern=character(0), fnFs, fReads)

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

plotQualityProfile(fnFs[1:25])
plotQualityProfile(fnFs[28:45])

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = FALSE)

# Have R get all primer orients 
# Create all orientations of the input sequence- we are only using the R1 (5' or FWD) primer here. 

FWD <- "CTTGGTCATTTAGAGGAAGTAA" 
REV <- "GCTGCGTTCTTCATCGATGC"

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = FALSE)

# Counts number of reads in which the primer is found
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]))

# NOTE- Miller human ITS data returned 145 seqs w/ fwd and 40 seqs w/ revcomp
# Skip cutadapt   

filtFs <- file.path(path, "09132021", basename(fnFs))
out <- filterAndTrim(fnFs, filtFs, truncQ = 2, maxEE = 5, truncLen = 150, 
                     rm.phix = TRUE, compress = TRUE, multithread = FALSE
)  # on windows, set multithread = FALSE
head(out)
plotQualityProfile(filtFs[1:18])
plotQualityProfile(filtFs[28:45])

# Check error rates and plot them

errF <- learnErrors(filtFs, multithread = FALSE)
plotErrors(errF, nominalQ = TRUE)

# Derep + DADA2 algorithm

derepFs <- derepFastq(filtFs, verbose = TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
dadaFs <- dada(derepFs, err = errF, multithread = FALSE)
dadaFs[[1]]       # Gives you ASV numbers from each sample, [1]= sample 1.

# Make sequencing table

seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# Remove Chimeras and make a new table

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
dim(seqtab.nochim)
table(nchar(getSequences(seqtab.nochim)))

# Track reads through the pipeline

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(fnFs, getN),sapply(filtFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("raw", "filt", "final")
rownames(track) <- sample.names
head(track)

# Seqtab.nochim is basically the "OTU" table. You can write the Fasta and "OTU" counts. 
# Give our seq headers more manageable names (HT_ASV_1, HT_ASV_2...)

asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")
for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">HTpanc_ASV", i, sep="_")
}


# Make the ASV count table:

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "HTpanc_counts.tsv", sep="\t", quote=F, col.names=NA)

######################### Assign Taxonomy #######################################

unite.ref <- "C:/Users/af206/Desktop/09082021_AF_DADA2/sh_general_release_dynamic_s_all_02.02.2019.fasta"  ## CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = FALSE, tryRC = TRUE, verbose = TRUE)
taxa.print <- taxa  ## Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

# Write a tax table:

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "HTpanc_tax.tsv", sep="\t", quote=F, col.names=NA)

# Making and writing out a fasta of our final ASV seqs:

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "09132021_HTpanc.fa")

