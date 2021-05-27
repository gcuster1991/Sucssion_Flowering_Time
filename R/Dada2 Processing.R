#dada2 processing
#Microbe Succession vs. Development Analyses
rm(list = ls())
setwd("~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/ATFT_merged_fastq")
#Acquire dada2, which is a bioconducter package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("xcode")

#Install Bioconductor packages 
BiocManager::install("dada2", version = "3.10")
BiocManager::install("S4Vectors")
BiocManager::install("DECIPHER")
BiocManager::install("phyloseq")
install.packages("rlang")

install.packages("devtools")

library(dada2)
library(phyloseq)
library("devtools")
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn)
library(Biostrings); packageVersion("Biostrings")
library(ggplot2); packageVersion("ggplot2")

#cp rdibner@teton.arcc.uwyo.edu:/project/brqtl/
# File parsing
#Notes from Gordon in Linda's class: DADA2 assumes the reads are split into forward and reverse. You will want to make two folders, one named "forward" and the other "reverse". You want to place all forward reads in the forward folder and the same for the reverse reads. The first two lines here tell R where it can find the reads. 
pathF <- "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/ATFT_merged_fastq/Forward"  # CHANGE ME to the directory containing your demultiplexed forward-read fastq files
pathR <- "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/ATFT_merged_fastq/Reverse"  # CHANGE ME TOO!
list.files(pathF)

#This next section createds a new folder inside the forward and reverse folders for the filtered reads. This is where you will be working from once you have filtered and trimmed the raw reads. 
filtpathF <- file.path(pathF, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...
fastqFs <- sort(list.files(pathF, pattern="R1.fq", full.names = T))
fastqRs <- sort(list.files(pathR, pattern="R2.fq", full.names = T))

quartz()
par(mfrow=c(1,2))
plotQualityProfile(fastqFs[1:6])
plotQualityProfile(fastqRs[1:6])

#You will have to rerun this in order to get the filter and trip parameter to work. This is an artifact of plotting the quality scores. 
fastqFs <- sort(list.files(pathF, pattern="R1.fq"))
fastqRs <- sort(list.files(pathR, pattern="R2.fq"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#From Gordon: This is the first step we have encountered in the DADA2 pipeline that has many options for customization and as such, the first step which requires you to fully understand what the options mean and to implement them. The truncLen parameter will cut your reads at a set length. The first number corresponds to the forward reads, and the second number corresponds to the reverse reads. In addition to cutting your reads at a set length you can also use a minimum quality score. To apply a minimum quality score, use truncQ. 

track.filt<-filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE, compress=TRUE, verbose=TRUE, multithread=TRUE)


# File parsing
filtpathF <- "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/ATFT_merged_fastq/Forward/filtered/" # CHANGE ME to the directory containing your filtered forward fastq files
filtpathR <- "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/ATFT_merged_fastq/Reverse/filtered/" # CHANGE ME ...
filtFs <- list.files(filtpathF, pattern="R1.fq", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="R2.fq", full.names = TRUE)
sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz. note to self: REILLY--KEEP FILES ZIPPED!!!
#The ! reads is 'NOT.' So this line below reads "if the sample.names and sample.namesR are not identical than stop"
#This is a final check to ensure the data exists in forward and reverse reads for the same samples. 
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names


#The learnErrors() step uses machine learning to uncover the true underlying proportion of errors. This is accomplished by alternating the estimation of the error rate and the inference of the sample composition until they converge.

#set.seed is used for reproducibilty 
set.seed(100)
# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)
#these plots allow you to see the error rates for each possible transition. 
plotErrors(errF)
plotErrors(errR)


#This next step merges forward and reverse sequences into a single read and keeps track of the abundance of each read. 

#This step also helps to reduce computing time; in that it removes duplicate sequences, greatly reducing the size of the data set. The dada() step infers the number of truly unique sequences present in the data set. This, along with the learnErrors() step, is the most important and unique step of the DADA2 pipeline.

#Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))

names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

#For the next step, we take the merged reads from above and move the data into a sequence table. This is the classic format for storing processed HTS data. A sequence variant table is compact and contains the information necessary for further analysis. No matter which pipeline you use, you will come to a step such as this. Other pipeliness call this table an OTU table because instead of sequence variants, they use OTUs grouped at some level of sequence similarity. The final product is the same, a table consisting of sequence counts by sites. 
seqtab <- makeSequenceTable(mergers)
#Remove chimeras
seqtabNoC <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)

#Track the number of reads retained after each step throughout the pipeline. This chunk of code below allows you to track the number of reads. There is no magic number or percentage of reads retained and you will want to address this on a data set by data set basis. 

#remember the Unique() function? This is another variation of it. 
getNreads <- function(x) sum(getUniques(x))
track <- cbind(track.filt, sapply(mergers, getNreads), rowSums(seqtab), rowSums(seqtabNoC))
colnames(track) <- c("input", "filtered", "merged", "tabled", "nonchim")
track
colSums(track)

#Now that we have our sequence table and have seen that each column is headed by the sequence that represents that variant, we want to assign taxonomy. In this step, we tell the assignTaxonomy() function several things. The first, which sequence table to use. We want it to use the table with no chimeras as this is our best "guess" at the true underlying proportions of microbes. Second, we tell the function where to find the database. Many different databases exist so you will want to figure out which one best suits your needs. Information on available databases can be found at https://benjjneb.github.io/dada2/training.html. 


# Assign taxonomy
taxa <- assignTaxonomy(seqtabNoC, "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/silva_nr_v132_train_set.fa.gz", minBoot = 70, multithread=TRUE)

#Optional additional step
taxa.species <- addSpecies(taxa, "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/silva_species_assignment_v132.fa.gz", tryRC=TRUE, verbose=TRUE)


#View taxa, and remove row names first:
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


saveRDS(seqtab, "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/seqtab_final_rrd.rds")
saveRDS(taxa, "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/tax_final_rrd.rds")
saveRDS(taxa.species, "~/Desktop/SuccessionVDevelopment/Flowering_Time_Sequence/tax_sp_final_rrd.rds")

#In the end, we have filtered and trimmed our data, learned errors, removed and or fixed non-biologically relevant variants, removed chimeras, created a sequence variant table, and assigned taxonomy. We now have all the necessary pieces for downstream analysis. The next tutorial utilizes the packages Phyloseq and Vegan for downstream analysis of our data. 

#theme_set(theme_bw())

#row.names(seqtabNoC)
#early.phenol <- which(row.names(seqtabNoC)==".early.")

