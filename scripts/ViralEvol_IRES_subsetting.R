# Read in packages
install.packages("Biostrings")
install.packages("phylotools")
install.packages("tidyverse")
library(Biostrings)
library(phylotools)
library(dplyr)
library(tidyr)
library(stringr)

# Set wd
setwd("~/Desktop/ViralAncestral/data")

# Use 5UTR' Region nucleotide data to filter individual sequences
fiveUTR <- readDNAStringSet("./5UTR_Region/picornaviridae2_5UTR_with_outgroup.fasta")
seq_name <- names(fiveUTR)
sequence <- paste(fiveUTR)
df <- data.frame(seq_name, sequence)

# Identify IRES subtypes
ires <- read.table("IRES_list.txt", fill = TRUE, sep = ",", header = FALSE)
ires <- ires %>% separate(V2, c("Name1", "Name2", "Name3", "Name4"))

# Manually impute values for obs that have excess information
ires$Name4[5] <- "IRES5"
ires$Name4[23] <- "IRES1"
ires$Name4[33] <- "IRES0"
ires$Name4[34] <- "IRES0"
ires$Name4[47] <- "IRES0"
ires$Name4[62] <- "IRES2"
ires$Name4[100] <- "IRES2"

# Filtering IRES subtypes
ires0 <- ires %>% filter(Name4 == "IRES0")
ires0$V1 <-gsub(">","",as.character(ires0$V1))

ires1 <- ires %>% filter(Name4 == "IRES1")
ires1$V1 <-gsub(">","",as.character(ires1$V1))

ires2 <- ires %>% filter(Name4 == "IRES2")
ires2$V1 <-gsub(">","",as.character(ires2$V1))

ires3 <- ires %>% filter(Name4 == "IRES3")
ires3$V1 <-gsub(">","",as.character(ires3$V1))

ires4 <- ires %>% filter(Name4 == "IRES4")
ires4$V1 <-gsub(">","",as.character(ires4$V1))

ires5 <- ires %>% filter(Name4 == "IRES5")
ires5$V1 <-gsub(">","",as.character(ires5$V1))

# Make data frames of IRES subtypes
df[c("Main", "Other")] <- str_split_fixed(df$seq_name, '_', 2)
df[c("Other", "Extra")] <- str_split_fixed(df$Other, ',', 2)

ires0[c('Main', 'Other')] <- str_split_fixed(ires0$V1, '_', 2)
df0 <- df %>% filter(Other %in% ires0$Other)
df0 <- df0[,c(1:2)]

ires1[c('Main', 'Other')] <- str_split_fixed(ires1$V1, '_', 2)
df1 <- df %>% filter(Other %in% ires1$Other)
df1 <- df1[,c(1:2)]

ires2[c('Main', 'Other')] <- str_split_fixed(ires2$V1, '_', 2)
df2 <- df %>% filter(Other %in% ires2$Other)
df2 <- df2[,c(1:2)]

ires3[c('Main', 'Other')] <- str_split_fixed(ires3$V1, '_', 2)
df3 <- df %>% filter(Other %in% ires3$Other)
df3 <- df3[,c(1:2)]

ires4[c('Main', 'Other')] <- str_split_fixed(ires4$V1, '_', 2)
df4 <- df %>% filter(Other %in% ires4$Other)
df4 <- df4[,c(1:2)]

ires5[c('Main', 'Other')] <- str_split_fixed(ires5$V1, '_', 2)
df5 <- df %>% filter(Other %in% ires5$Other)
df5 <- df5[,c(1:2)]

# Make new df to create fasta file outputs
newdf <- read.fasta("./5UTR_Region/picornaviridae2_5UTR_with_outgroup.fasta")

# fasta files for IRES type from 5'UTR region
# location of sequences the same between Replicase Gene and 5'UTR Region - so will use for replicase gene too with respect to filtering 
ires_0 <- newdf %>% filter(seq_name %in% df0$seq_name)
dat2fasta(ires_0, outfile = "5UTR_ires0_NT.fasta")

ires_1 <- newdf %>% filter(seq_name %in% df1$seq_name)
dat2fasta(ires_1, outfile = "5UTR_ires1_NT.fasta")

ires_2 <- newdf %>% filter(seq_name %in% df2$seq_name)
dat2fasta(ires_2, outfile = "5UTR_ires2_NT.fasta")

ires_3 <- newdf %>% filter(seq_name %in% df3$seq_name)
dat2fasta(ires_3, outfile = "5UTR_ires3_NT.fasta")

ires_4 <- newdf %>% filter(seq_name %in% df4$seq_name)
dat2fasta(ires_4, outfile = "5UTR_ires4_NT.fasta")

ires_5 <- newdf %>% filter(seq_name %in% df5$seq_name)
dat2fasta(ires_5, outfile = "5UTR_ires5_NT.fasta")

## Read in Nucleotides for Replicase Gene
RepNT <- read.fasta("./Replicase_Gene/PRNT.fasta")
RepAA <- read.fasta("./Replicase_Gene/PRAA.fasta")
RepNT_copy <- RepNT
RepAA_copy <- RepAA

## Create data frames for each IRES type - manual identification from previous IRES subsetting
Replicase_ires0_AA <- RepAA_copy[c(6, 7, 8, 22, 23,35, 36, 37, 38, 40, 42, 47, 51, 61, 62, 65, 66, 68, 69),]
Replicase_ires0_NT <- RepNT_copy[c(6, 7, 8, 22, 23,35, 36, 37, 38, 40, 42, 47, 51, 61, 62, 65, 66, 68, 69),]

RepAA_copy <- RepAA_copy[-c(6, 7, 8, 22, 23,35, 36, 37, 38, 40, 42, 47, 51, 61, 62, 65, 66, 68, 69),]
RepNT_copy <- RepNT_copy[-c(6, 7, 8, 22, 23,35, 36, 37, 38, 40, 42, 47, 51, 61, 62, 65, 66, 68, 69),]

row.names(RepAA_copy) <- 1:nrow(RepAA_copy)
row.names(RepNT_copy) <- 1:nrow(RepNT_copy)

ires1_for_repAA <- RepAA_copy[c(1, 2, 4, 5, 6, 7, 8, 10, 27, 38, 43, 47, 52, 55, 64, 65),]
ires1_for_repNT <- RepNT_copy[c(1, 2, 4, 5, 6, 7, 8, 10, 27, 38, 43, 47, 52, 55, 64, 65),]

RepAA_copy <- RepAA_copy[-c(1, 2, 4, 5, 6, 7, 8, 10, 27, 38, 43, 47, 52, 55, 64, 65),]
RepNT_copy <- RepNT_copy[-c(1, 2, 4, 5, 6, 7, 8, 10, 27, 38, 43, 47, 52, 55, 64, 65),]

row.names(RepAA_copy) <- 1:nrow(RepAA_copy)
row.names(RepNT_copy) <- 1:nrow(RepNT_copy)

ires2_for_repAA <- RepAA_copy[c(1, 3, 6, 8, 11, 12, 13, 14, 15, 18, 22, 23, 26, 28, 32, 36, 37, 47),]
ires2_for_repNT <- RepNT_copy[c(1, 3, 6, 8, 11, 12, 13, 14, 15, 18, 22, 23, 26, 28, 32, 36, 37, 47),]

RepAA_copy <- RepAA_copy[-c(1, 3, 6, 8, 11, 12, 13, 14, 15, 18, 22, 23, 26, 28, 32, 36, 37, 47),]
RepNT_copy <- RepNT_copy[-c(1, 3, 6, 8, 11, 12, 13, 14, 15, 18, 22, 23, 26, 28, 32, 36, 37, 47),]

row.names(RepAA_copy) <- 1:nrow(RepAA_copy)
row.names(RepNT_copy) <- 1:nrow(RepNT_copy)

ires3_for_repAA <- RepAA_copy[c(8, 25, 26),]
ires3_for_repNT <- RepNT_copy[c(8, 25, 26),]

RepAA_copy <- RepAA_copy[-c(8, 25, 26),]
RepNT_copy <- RepNT_copy[-c(8, 25, 26),]

row.names(RepAA_copy) <- 1:nrow(RepAA_copy)
row.names(RepNT_copy) <- 1:nrow(RepNT_copy)

ires4_for_repAA <- RepAA_copy[c(4, 5, 7, 10, 11, 12, 16, 20, 21),]
ires4_for_repNT <- RepNT_copy[c(4, 5, 7, 10, 11, 12, 16, 20, 21),]

RepAA_copy <- RepAA_copy[-c(4, 5, 7, 10, 11, 12, 16, 20, 21),]
RepNT_copy <- RepNT_copy[-c(4, 5, 7, 10, 11, 12, 16, 20, 21),]

row.names(RepAA_copy) <- 1:nrow(RepAA_copy)
row.names(RepNT_copy) <- 1:nrow(RepNT_copy)

ires5_for_repAA <- RepAA_copy[c(1, 7, 9),]
ires5_for_repNT <- RepNT_copy[c(1, 7, 9),]

RepAA_copy <- RepAA_copy[-c(1, 7, 9),]
RepNT_copy <- RepNT_copy[-c(1, 7, 9),]

row.names(RepAA_copy) <- 1:nrow(RepAA_copy)
row.names(RepNT_copy) <- 1:nrow(RepNT_copy)

# Make fasta files for each well-studied IRES type
dat2fasta(ires1_for_repAA, outfile = "Replicase_ires1_AA.fasta")
dat2fasta(ires1_for_repNT, outfile = "Replicase_ires1_NT.fasta")
dat2fasta(ires2_for_repAA, outfile = "Replicase_ires2_AA.fasta")
dat2fasta(ires2_for_repNT, outfile = "Replicase_ires2_NT.fasta")
dat2fasta(ires3_for_repAA, outfile = "Replicase_ires3_AA.fasta")
dat2fasta(ires3_for_repNT, outfile = "Replicase_ires3_NT.fasta")




