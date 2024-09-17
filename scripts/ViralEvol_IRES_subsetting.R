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
setwd("~/Desktop/ViralAncestral/Results_Fall_2024")

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
RepAA[85,3] <- "NC_004046"
RepAA[86,3] <- "NC_004052"
ires$id <- gsub(".*(NC_\\w+).*", "\\1", ires$V1)

# Focus on IRES 1 
# 5'UTR sequences
ires1_5UTR <- fiveUTR %>% filter(id %in% ires[ires$Name4=="IRES1",]$id) %>% select(-'id') # get only IRES 1
ires1_5UTR_with_outgroup <- rbind(ires1_5UTR,fiveUTR[130:131,-3]) # attach mitovirus outgroup 
ires1_5UTR_with_outgroup$seq.name <- sub("^[^ ]+ ", "", ires1_5UTR_with_outgroup$seq.name) # remove sequence code
ires1_5UTR_with_outgroup$seq.name <- gsub("complete genome", "", ires1_5UTR_with_outgroup$seq.name) # remove 'complete genome'
ires1_5UTR_with_outgroup$seq.name<- gsub(" ", "_", ires1_5UTR_with_outgroup$seq.name) # add underscores to retain in fasta file
dat2fasta(ires1_5UTR_with_outgroup,"~/Desktop/ViralAncestral/Results_Fall_2024/ires1_5UTR_with_outgroup.fasta")

# Replicase NT
ires1_Rep_NT <- RepNT %>% filter(id %in% ires[ires$Name4=="IRES1",]$id) %>% select(-'id')
ires1_Rep_NT_with_outgroup <- rbind(ires1_Rep_NT,RepNT[85:86,-3])
ires1_Rep_NT_with_outgroup$seq.name <- sub("^[^ ]+ ", "", ires1_Rep_NT_with_outgroup$seq.name) # remove sequence code
ires1_Rep_NT_with_outgroup$seq.name <- gsub("complete genome", "", ires1_Rep_NT_with_outgroup$seq.name) # remove 'complete genome'
ires1_Rep_NT_with_outgroup$seq.name<- gsub(" ", "_", ires1_Rep_NT_with_outgroup$seq.name) # add underscores to retain in fasta file
dat2fasta(ires1_Rep_NT_with_outgroup ,"~/Desktop/ViralAncestral/Results_Fall_2024/ires1_REP_NT_with_outgroup.fasta")

# Replicase AA 
ires1_Rep_AA <- RepAA %>% filter(id %in% ires[ires$Name4=="IRES1",]$id) %>% select(-'id')
ires1_Rep_AA_with_outgroup <- rbind(ires1_Rep_AA,RepAA[85:86,-3])
ires1_Rep_AA_with_outgroup$seq.name <- sub("^[^ ]+ [^ ]+ [^ ]+ ", "", ires1_Rep_AA_with_outgroup$seq.name) # remove sequence code
ires1_Rep_AA_with_outgroup[17,]$seq.name <- ires1_Rep_NT_with_outgroup[17,]$seq.name
ires1_Rep_AA_with_outgroup[18,]$seq.name <- ires1_Rep_NT_with_outgroup[18,]$seq.name
ires1_Rep_AA_with_outgroup$seq.name <- gsub("complete genome", "", ires1_Rep_AA_with_outgroup$seq.name) # remove 'complete genome'
ires1_Rep_AA_with_outgroup$seq.name<- gsub(" ", "_", ires1_Rep_AA_with_outgroup$seq.name)
dat2fasta(ires1_Rep_AA_with_outgroup ,"~/Desktop/ViralAncestral/Results_Fall_2024/ires1_REP_AA_with_outgroup.fasta")



## Iterate over IRES types from 0 to 5 to filter by IRES type and write to FASTA file ## 

# 5'UTR Region 
for (i in 1:1) { ## only working with IRES 1 right now to get df to attach mitovirus
  # Filter data based on IRES type
  filtered_df <- fiveUTR %>%
    filter(id %in% ires$id[ires$Name4 == paste0("IRES", i)]) #%>% # filter by IRES type i 
    #select(seq.name, seq.text) # select name and sequence --> need these column names for dat2fasta
  
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

# Add mitoviruses to IRES 1
filtered_df[23:24,] <- fiveUTR[130:131, 1:2]
# Separate name from sequence
filtered_df[c('id', 'name')] <- str_split_fixed(filtered_df$seq.name, ' ', 2)
# put underscores so they're uniquely identified in downstream analyses
filtered_df$name <- gsub(" ", "_", filtered_df$name)

# Clean data so it's IRES 1 virus names with sequences
IRES1_with_outgroup <- data.frame(cbind(filtered_df$name, filtered_df$seq.text))
IRES1_with_outgroup[c('name', 'type')] <- str_split_fixed(IRES1_with_outgroup$X1, ',', 2)
IRES1_with_outgroup <- data.frame(cbind(IRES1_with_outgroup$name, IRES1_with_outgroup$X2))
colnames(IRES1_with_outgroup) <- c("seq.name", "seq.text")
new_outfile_path <- sprintf("~/Desktop/ViralAncestral/data/5UTR_Region/IRES_Subtypes/IRES1_with_outgroup.fasta")
dat2fasta(IRES1_with_outgroup, outfile = new_outfile_path)

# Rearrange columns and remove original name column
df <- df[c('First Name', 'Last Name', 'State')]


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
