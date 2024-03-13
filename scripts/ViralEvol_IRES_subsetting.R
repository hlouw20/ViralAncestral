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
setwd("~/Desktop/ViralAncestral/data/Original Data March 2023 (Start Here)")

# Read in 5'UTR NT, Replicase NT, Replicase AA, and IRES types
fiveUTR <- read.fasta("picornaviridae_5UTR_with_outgroup.fasta")
RepNT <- read.fasta("picornaviridae_replicase_nucleotide_with_outgroup.fasta")
RepAA <- read.fasta("picornaviridae_replicase_aminoacid_with_outgroup.fasta")
ires <- read.table("IRES_list.txt", fill = TRUE, sep = ",", header = FALSE)

# Isolate IRES types and manually impute missing values 
ires <- ires %>% separate(V2, c("Name1", "Name2", "Name3", "Name4"))
ires$Name4[5] <- "IRES5"
ires$Name4[23] <- "IRES1"
ires$Name4[33] <- "IRES0"
ires$Name4[34] <- "IRES0"
ires$Name4[47] <- "IRES0"
ires$Name4[62] <- "IRES2"
ires$Name4[100] <- "IRES2"

# Extract sequence IDs for filtering
fiveUTR$id <- gsub(".*(NC_\\w+).*", "\\1", fiveUTR$seq.name)
RepNT$id <- gsub(".*(NC_\\w+).*", "\\1", RepNT$seq.name)
RepAA$id <- gsub(".*(NC_\\w+).*", "\\1", RepAA$seq.name)
ires$id <- gsub(".*(NC_\\w+).*", "\\1", ires$V1)

## Iterate over IRES types from 0 to 5 to filter by IRES type and write to FASTA file ## 

# 5'UTR Region 
for (i in 0:5) {
  # Filter data based on IRES type
  filtered_df <- fiveUTR %>%
    filter(id %in% ires$id[ires$Name4 == paste0("IRES", i)]) %>% # filter by IRES type i 
    select(seq.name, seq.text) # select name and sequence --> need these column names for dat2fasta
  
  # Check if filtered_df is empty
  if (nrow(filtered_df) == 0) {
    cat("No data found for IRES", i, "\n")
  } else {
    # Define the output file path
    outfile_path <- sprintf("~/Desktop/ViralAncestral/data/5UTR_Region/IRES_Subtypes/5UTR_ires%d_NT.fasta", i)
    
    # Call dat2fasta with filtered data and outfile path
    dat2fasta(filtered_df, outfile = outfile_path) # write to FASTA file 
  }
}


# Replicase NT
for (i in 0:5) {
  # Filter data based on IRES type
  filtered_df <- RepNT %>%
    filter(id %in% ires$id[ires$Name4 == paste0("IRES", i)]) %>% # filter by IRES type i 
    select(seq.name, seq.text)  # select name and sequence --> need these column names for dat2fasta
  
  # Check if filtered_df is empty
  if (nrow(filtered_df) == 0) {
    cat("No data found for IRES", i, "\n")
  } else {
    # Define the output file path
    outfile_path <- sprintf("~/Desktop/ViralAncestral/data/Replicase_Gene/IRES_Subtypes/Replicase_ires%d_NT.fasta", i)
    
    # Call dat2fasta with filtered data and outfile path
    dat2fasta(filtered_df, outfile = outfile_path) # write to FASTA file
  }
}

# Replicase AA
for (i in 0:5) {
  # Filter data based on IRES type
  filtered_df <- RepAA %>%
    filter(id %in% ires$id[ires$Name4 == paste0("IRES", i)]) %>% # filter by IRES type i 
    select(seq.name, seq.text) # select name and sequence --> need these column names for dat2fasta
  
  # Check if filtered_df is empty
  if (nrow(filtered_df) == 0) {
    cat("No data found for IRES", i, "\n")
  } else {
    # Define the output file path
    outfile_path <- sprintf("~/Desktop/ViralAncestral/data/Replicase_Gene/IRES_Subtypes/Replicase_ires%d_AA.fasta", i)
    
    # Call dat2fasta with filtered data and outfile path
    dat2fasta(filtered_df, outfile = outfile_path) # write to FASTA file
  }
}
