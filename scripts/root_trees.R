## Root Trees from IQ-Tree using Midpoint Rooting ## 

library(phytools)
setwd("~/Desktop/ViralAncestral/results")

## 5'UTR Trees 

## MUSCLE only alignment under GTR+F+I+G4
muscle_gtr_5utr_nt_ires1 <- read.tree("./5UTR_Trees/GTR+F+I+G4/ML_Tree/5UTR_MUSCLE_ires1_GTR_NT.fasta.treefile")
muscle_gtr_5utr_nt_ires2 <- read.tree("./5UTR_Trees/GTR+F+I+G4/ML_Tree/5UTR_MUSCLE_ires2_GTR_NT.fasta.treefile")
muscle_gtr_5utr_nt_ires3 <- read.tree("./5UTR_Trees/GTR+F+I+G4/ML_Tree/5UTR_MUSCLE_ires3_GTR_NT.fasta.treefile")

muscle_gtr_5utr_nt_ires1_rooted <- midpoint.root(muscle_gtr_5utr_nt_ires1, node.lables = "support")
muscle_gtr_5utr_nt_ires2_rooted <- midpoint.root(muscle_gtr_5utr_nt_ires2, node.lables = "support")
muscle_gtr_5utr_nt_ires3_rooted <- midpoint.root(muscle_gtr_5utr_nt_ires3, node.lables = "support")

write.tree(muscle_gtr_5utr_nt_ires1_rooted, file="./Rooted_Trees/5UTR/muscle_gtr_5utr_nt_ires1_rooted.treefile")
write.tree(muscle_gtr_5utr_nt_ires2_rooted, file="./Rooted_Trees/5UTR/muscle_gtr_5utr_nt_ires2_rooted.treefile")
write.tree(muscle_gtr_5utr_nt_ires3_rooted, file="./Rooted_Trees/5UTR/muscle_gtr_5utr_nt_ires3_rooted.treefile")

## MUSCLE only alignment under Extended Model Finder
muscle_sym_5utr_nt_ires1 <- read.tree("./5UTR_Trees/MFP/ML_Tree/5UTR_MUSCLE_ires1_MFP_NT.fasta.treefile")
muscle_tpm3u_5utr_nt_ires2 <- read.tree("./5UTR_Trees/MFP/ML_Tree/5UTR_MUSCLE_ires2_MFP_NT.fasta.treefile")
muscle_tn_5utr_nt_ires3 <- read.tree("./5UTR_Trees/MFP/ML_Tree/5UTR_MUSCLE_ires3_MFP_NT.fasta.treefile")

muscle_sym_5utr_nt_ires1_rooted  <- midpoint.root(muscle_sym_5utr_nt_ires1, node.lables = "support")
muscle_tpm3u_5utr_nt_ires2_rooted <- midpoint.root(muscle_tpm3u_5utr_nt_ires2, node.lables = "support")
muscle_tn_5utr_nt_ires3_rooted <- midpoint.root(muscle_tn_5utr_nt_ires3, node.labels = "support")

write.tree(muscle_sym_5utr_nt_ires1_rooted, file="./Rooted_Trees/5UTR/muscle_sym_5utr_nt_ires1_rooted.treefile")
write.tree(muscle_tpm3u_5utr_nt_ires2_rooted, file="./Rooted_Trees/5UTR/muscle_tpm3u_5utr_nt_ires2_rooted.treefile")
write.tree(muscle_tn_5utr_nt_ires3_rooted, file="./Rooted_Trees/5UTR/muscle_tn_5utr_nt_ires3_rooted.treefile")

## MUSCLE + HmmCleaner alignment under GTR+F+I+G4 
muscle_hmm_gtr_5utr_nt_ires1 <- read.tree("./5UTR_Trees/GTR+F+I+G4/ML_Tree/5UTR_MUSCLE_ires1_GTR_NT_hmm.fasta.treefile")
muscle_hmm_gtr_5utr_nt_ires2 <- read.tree("./5UTR_Trees/GTR+F+I+G4/ML_Tree/5UTR_MUSCLE_ires2_GTR_NT_hmm.fasta.treefile")
muscle_hmm_gtr_5utr_nt_ires3 <- read.tree("./5UTR_Trees/GTR+F+I+G4/ML_Tree/5UTR_MUSCLE_ires3_GTR_NT_hmm.fasta.treefile")

muscle_hmm_gtr_5utr_nt_ires1_rooted <- midpoint.root(muscle_hmm_gtr_5utr_nt_ires1, node.lables = "support")
muscle_hmm_gtr_5utr_nt_ires2_rooted <- midpoint.root(muscle_hmm_gtr_5utr_nt_ires2, node.lables = "support")
muscle_hmm_gtr_5utr_nt_ires3_rooted <- midpoint.root(muscle_hmm_gtr_5utr_nt_ires3, node.lables = "support")

write.tree(muscle_hmm_gtr_5utr_nt_ires1_rooted, file="./Rooted_Trees/5UTR/muscle_hmm_gtr_5utr_nt_ires1_rooted.treefile")
write.tree(muscle_hmm_gtr_5utr_nt_ires2_rooted, file="./Rooted_Trees/5UTR/muscle_hmm_gtr_5utr_nt_ires2_rooted.treefile")
write.tree(muscle_hmm_gtr_5utr_nt_ires3_rooted, file="./Rooted_Trees/5UTR/muscle_hmm_gtr_5utr_nt_ires3_rooted.treefile")

## MUSCLE + HmmCleaner alignment under Extended Model Finder 
muscle_hmm_sym_5utr_nt_ires1 <- read.tree("./5UTR_Trees/MFP/ML_Tree/5UTR_MUSCLE_ires1_MFP_NT_hmm.fasta.treefile")
muscle_hmm_tim3e_5utr_nt_ires2 <- read.tree("./5UTR_Trees/MFP/ML_Tree/5UTR_MUSCLE_ires2_MFP_NT_hmm.fasta.treefile")
muscle_hmm_hky_5utr_nt_ires3 <- read.tree("./5UTR_Trees/MFP/ML_Tree/5UTR_MUSCLE_ires3_MFP_NT_hmm.fasta.treefile")

muscle_hmm_sym_5utr_nt_ires1_rooted  <- midpoint.root(muscle_hmm_sym_5utr_nt_ires1, node.lables = "support")
muscle_hmm_tim3e_5utr_nt_ires2_rooted <- midpoint.root(muscle_hmm_tim3e_5utr_nt_ires2, node.lables = "support")
muscle_hmm_hky_5utr_nt_ires3_rooted <- midpoint.root(muscle_hmm_hky_5utr_nt_ires3, node.labels = "support")

write.tree(muscle_hmm_sym_5utr_nt_ires1_rooted, file="./Rooted_Trees/5UTR/muscle_hmm_sym_5utr_nt_ires1_rooted.treefile")
write.tree(muscle_hmm_tim3e_5utr_nt_ires2_rooted, file="./Rooted_Trees/5UTR/muscle_hmm_tim3e_5utr_nt_ires2_rooted.treefile")
write.tree(muscle_hmm_hky_5utr_nt_ires3_rooted, file="./Rooted_Trees/5UTR/muscle_hmm_hky_5utr_nt_ires3_rooted.treefile")

## Replicase Trees ##

## MUSCLE NT only alignment under GTR+F+I+G4
muscle_gtr_rep_nt_ires1 <- read.tree("./Replicase_Trees/GTR+F+I+G4/ML_Tree/Rep_MUSCLE_ires1_GTR_NT.fasta.treefile")
muscle_gtr_rep_nt_ires2 <- read.tree("./Replicase_Trees/GTR+F+I+G4/ML_Tree/Rep_MUSCLE_ires2_GTR_NT.fasta.treefile")
muscle_gtr_rep_nt_ires3 <- read.tree("./Replicase_Trees/GTR+F+I+G4/ML_Tree/Rep_MUSCLE_ires3_GTR_NT.fasta.treefile")

muscle_gtr_rep_nt_ires1_rooted <- midpoint.root(muscle_gtr_rep_nt_ires1, node.lables = "support")
muscle_gtr_rep_nt_ires2_rooted <- midpoint.root(muscle_gtr_rep_nt_ires2, node.lables = "support")
muscle_gtr_rep_nt_ires3_rooted <- midpoint.root(muscle_gtr_rep_nt_ires3, node.lables = "support")

write.tree(muscle_gtr_rep_nt_ires1_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_gtr_rep_nt_ires1_rooted.treefile")
write.tree(muscle_gtr_rep_nt_ires2_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_gtr_rep_nt_ires2_rooted.treefile")
write.tree(muscle_gtr_rep_nt_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_gtr_rep_nt_ires3_rooted.treefile")

## MUSCLE NT only alignment under Extended Model Finder --> only need to do for IRES 3
muscle_tpm2u_rep_nt_ires3 <- read.tree("./Replicase_Trees/MFP/ML_Tree/Rep_MUSCLE_ires3_MFP_NT.fasta.treefile")
muscle_tpm2u_rep_nt_ires3_rooted <- midpoint.root(muscle_tpm2u_rep_nt_ires3, node.lables = "support")
write.tree(muscle_tpm2u_rep_nt_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_tpm2u_rep_nt_ires3_rooted.treefile")

## MUSCLE AA only alignment under LG+F+I+G4
muscle_lg_rep_aa_ires1 <- read.tree("./Replicase_Trees/LG+F+I+G4/ML_Tree/Rep_MUSCLE_ires1_LG_AA.fasta.treefile")
muscle_lg_rep_aa_ires2 <- read.tree("./Replicase_Trees/LG+F+I+G4/ML_Tree/Rep_MUSCLE_ires2_LG_AA.fasta.treefile")
muscle_lg_rep_aa_ires3 <- read.tree("./Replicase_Trees/LG+F+I+G4/ML_Tree/Rep_MUSCLE_ires3_LG_AA.fasta.treefile")

muscle_lg_rep_aa_ires1_rooted <- midpoint.root(muscle_lg_rep_aa_ires1, node.lables = "support")
muscle_lg_rep_aa_ires2_rooted <- midpoint.root(muscle_lg_rep_aa_ires2, node.lables = "support")
muscle_lg_rep_aa_ires3_rooted <- midpoint.root(muscle_lg_rep_aa_ires3, node.lables = "support")

write.tree(muscle_lg_rep_aa_ires1_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_lg_rep_aa_ires1_rooted.treefile")
write.tree(muscle_lg_rep_aa_ires2_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_lg_rep_aa_ires2_rooted.treefile")
write.tree(muscle_lg_rep_aa_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_lg_rep_aa_ires3_rooted.treefile")

## MUSCLE AA only alignment under FLAVI+F+I+G4
muscle_flavi_rep_aa_ires1 <- read.tree("./Replicase_Trees/FLAVI+F+I+G4/ML_Tree/Rep_MUSCLE_ires1_FLAVI_AA.fasta.treefile")
muscle_flavi_rep_aa_ires2 <- read.tree("./Replicase_Trees/FLAVI+F+I+G4/ML_Tree/Rep_MUSCLE_ires2_FLAVI_AA.fasta.treefile")
muscle_flavi_rep_aa_ires3 <- read.tree("./Replicase_Trees/FLAVI+F+I+G4/ML_Tree/Rep_MUSCLE_ires3_FLAVI_AA.fasta.treefile")

muscle_flavi_rep_aa_ires1_rooted <- midpoint.root(muscle_flavi_rep_aa_ires1, node.lables = "support")
muscle_flavi_rep_aa_ires2_rooted <- midpoint.root(muscle_flavi_rep_aa_ires2, node.lables = "support")
muscle_flavi_rep_aa_ires3_rooted <- midpoint.root(muscle_flavi_rep_aa_ires3, node.lables = "support")

write.tree(muscle_flavi_rep_aa_ires1_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_flavi_rep_aa_ires1_rooted.treefile")
write.tree(muscle_flavi_rep_aa_ires2_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_flavi_rep_aa_ires2_rooted.treefile")
write.tree(muscle_flavi_rep_aa_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_flavi_rep_aa_ires3_rooted.treefile")

## MUSCLE AA only alignment under Extended Model Finder --> only need to do for IRES 3
muscle_blosum62_rep_aa_ires3 <- read.tree("./Replicase_Trees/MFP/ML_Tree/Rep_MUSCLE_ires3_MFP_AA.fasta.treefile")
muscle_blosum62_rep_aa_ires3_rooted <- midpoint.root(muscle_blosum62_rep_aa_ires3, node.lables = "support")
write.tree(muscle_blosum62_rep_aa_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_blosum62_rep_aa_ires3_rooted.treefile")

## MUSCLE NT + HmmCleaner alignment under GTR+F+I+G4 Model 
muscle_hmm_gtr_rep_nt_ires1 <- read.tree("./Replicase_Trees/GTR+F+I+G4/ML_Tree/Rep_MUSCLE_ires1_GTR_NT_hmm.fasta.treefile")
muscle_hmm_gtr_rep_nt_ires2 <- read.tree("./Replicase_Trees/GTR+F+I+G4/ML_Tree/Rep_MUSCLE_ires2_GTR_NT_hmm.fasta.treefile")
muscle_hmm_gtr_rep_nt_ires3 <- read.tree("./Replicase_Trees/GTR+F+I+G4/ML_Tree/Rep_MUSCLE_ires3_GTR_NT_hmm.fasta.treefile")

muscle_hmm_gtr_rep_nt_ires1_rooted <- midpoint.root(muscle_hmm_gtr_rep_nt_ires1, node.lables = "support")
muscle_hmm_gtr_rep_nt_ires2_rooted <- midpoint.root(muscle_hmm_gtr_rep_nt_ires2, node.lables = "support")
muscle_hmm_gtr_rep_nt_ires3_rooted <- midpoint.root(muscle_hmm_gtr_rep_nt_ires3, node.lables = "support")

write.tree(muscle_hmm_gtr_rep_nt_ires1_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_gtr_rep_nt_ires1_rooted.treefile")
write.tree(muscle_hmm_gtr_rep_nt_ires2_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_gtr_rep_nt_ires2_rooted.treefile")
write.tree(muscle_hmm_gtr_rep_nt_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_gtr_rep_nt_ires3_rooted.treefile")

## MUSCLE NT + HmmCleaner alignment under Extended Model Finder --> Only need to do IRES 3
muscle_hmm_tpm2u_rep_nt_ires3 <- read.tree("./Replicase_Trees/MFP/ML_Tree/Rep_MUSCLE_ires3_MFP_NT_hmm.fasta.treefile")
muscle_hmm_tpm2u_rep_nt_ires3_rooted <- midpoint.root(muscle_hmm_tpm2u_rep_nt_ires3, node.lables = "support")
write.tree(muscle_hmm_tpm2u_rep_nt_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_tpm2u_rep_nt_ires3_rooted.treefile")

## MUSCLE AA + HmmCleaner alignment under LG+F+I+G4 Model
muscle_hmm_lg_rep_aa_ires1 <- read.tree("./Replicase_Trees/LG+F+I+G4/ML_Tree/Rep_MUSCLE_ires1_LG_AA_hmm.fasta.treefile")
muscle_hmm_lg_rep_aa_ires2 <- read.tree("./Replicase_Trees/LG+F+I+G4/ML_Tree/Rep_MUSCLE_ires2_LG_AA_hmm.fasta.treefile")
muscle_hmm_lg_rep_aa_ires3 <- read.tree("./Replicase_Trees/LG+F+I+G4/ML_Tree/Rep_MUSCLE_ires3_LG_AA_hmm.fasta.treefile")

muscle_hmm_lg_rep_aa_ires1_rooted <- midpoint.root(muscle_hmm_lg_rep_aa_ires1, node.lables = "support")
muscle_hmm_lg_rep_aa_ires2_rooted <- midpoint.root(muscle_hmm_lg_rep_aa_ires2, node.lables = "support")
muscle_hmm_lg_rep_aa_ires3_rooted <- midpoint.root(muscle_hmm_lg_rep_aa_ires3, node.lables = "support")

write.tree(muscle_hmm_lg_rep_aa_ires1_rooted , file="./Rooted_Trees/Replicase_Gene/muscle_hmm_lg_rep_aa_ires1_rooted .treefile")
write.tree(muscle_hmm_lg_rep_aa_ires2_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_lg_rep_aa_ires2_rooted.treefile")
write.tree(muscle_hmm_lg_rep_aa_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_lg_rep_aa_ires3_rooted.treefile")

## MUSCLE AA + HmmCleaner alignment under FLAVI+F+I+G4 Model
muscle_hmm_flavi_rep_aa_ires1 <- read.tree("./Replicase_Trees/FLAVI+F+I+G4/ML_Tree/Rep_MUSCLE_ires1_FLAVI_AA_hmm.fasta.treefile")
muscle_hmm_flavi_rep_aa_ires2 <- read.tree("./Replicase_Trees/FLAVI+F+I+G4/ML_Tree/Rep_MUSCLE_ires2_FLAVI_AA_hmm.fasta.treefile")
muscle_hmm_flavi_rep_aa_ires3 <- read.tree("./Replicase_Trees/FLAVI+F+I+G4/ML_Tree/Rep_MUSCLE_ires3_FLAVI_AA_hmm.fasta.treefile")

muscle_hmm_flavi_rep_aa_ires1_rooted <- midpoint.root(muscle_hmm_flavi_rep_aa_ires1, node.lables = "support")
muscle_hmm_flavi_rep_aa_ires2_rooted <- midpoint.root(muscle_hmm_flavi_rep_aa_ires2, node.lables = "support")
muscle_hmm_flavi_rep_aa_ires3_rooted <- midpoint.root(muscle_hmm_flavi_rep_aa_ires3, node.lables = "support")

write.tree(muscle_hmm_flavi_rep_aa_ires1_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_flavi_rep_aa_ires1_rooted.treefile")
write.tree(muscle_hmm_flavi_rep_aa_ires2_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_flavi_rep_aa_ires2_rooted.treefile")
write.tree(muscle_hmm_flavi_rep_aa_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_flavi_rep_aa_ires3_rooted.treefile")

## MUSCLE AA + HmmCleaner alignment under Extended Model Finder --> only need to do IRES 3
muscle_hmm_blosum62_rep_aa_ires3 <- read.tree("./Replicase_Trees/MFP/ML_Tree/Rep_MUSCLE_ires3_MFP_AA_hmm.fasta.treefile")
muscle_hmm_blosum62_rep_aa_ires3_rooted <- midpoint.root(muscle_hmm_blosum62_rep_aa_ires3, node.lables = "support")
write.tree(muscle_hmm_blosum62_rep_aa_ires3_rooted, file="./Rooted_Trees/Replicase_Gene/muscle_hmm_blosum62_rep_aa_ires3_rooted.treefile")

