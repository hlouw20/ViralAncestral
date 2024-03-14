## Visualize Ancestral Trees ##

library(phytools)

## IRES 1 GTR - Untrimmed ## 

ires1_gtr_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Ancestral/GTR+I+G4/5UTR_MUSCLE_ires1_NT.fasta.raxml.ancestralTree")
ires1_gtr_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Ancestral/GTR+I+G4/fiveUTR_muscle_ires1_nt_gtr_ancestral.fasta")
plotTree(ires1_gtr_tree, ftype = "i", fize = 0.5, lwd = 1)
labelnodes(1:ires1_gtr_tree$Nnode+Ntip(ires1_gtr_tree),
           1:ires1_gtr_tree$Nnode+Ntip(ires1_gtr_tree),
           interactive = FALSE, cex=0.3, shape = "ellipse")

## IRES 1 GTR - Trimmed ## 

ires1_trimmed_gtr_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Trimmed_Ancestral/GTR+I+G4/5UTR_MUSCLE_ires1_NT_hmm.fasta.raxml.ancestralTree")
ires1_trimmed_gtr_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Trimmed_Ancestral/GTR+I+G4/fiveUTR_muscle_ires1_nt_hmm_gtr_ancestral.fasta")
plotTree(ires1_trimmed_gtr_tree, ftype = "i", fize = 0.5, lwd = 1)
labelnodes(1:ires1_trimmed_gtr_tree$Nnode+Ntip(ires1_trimmed_gtr_tree),
           1:ires1_trimmed_gtr_tree$Nnode+Ntip(ires1_trimmed_gtr_tree),
           interactive = FALSE, cex=0.3, shape = "ellipse")

## IRES 1 SYM - Untrimmed ## 

ires1_sym_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Ancestral/SYM+G4/5UTR_MUSCLE_ires1_NT.fasta.raxml.ancestralTree")
ires1_sym_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Ancestral/SYM+G4/fiveUTR_muscle_ires1_nt_sym_ancestral.fasta")
plotTree(ires1_sym_tree, ftype = "i", fize = 0.5, lwd = 1)
labelnodes(1:ires1_sym_tree$Nnode+Ntip(ires1_sym_tree),
           1:ires1_sym_tree$Nnode+Ntip(ires1_sym_tree),
           interactive = FALSE, cex=0.3, shape = "ellipse")

## IRES 1 SYM - Trimmed ## 

ires1_trimmed_sym_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Trimmed_Ancestral/SYM+G4/5UTR_MUSCLE_ires1_NT_hmm.fasta.raxml.ancestralTree")
ires1_trimmed_sym_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES1_MUSCLE_Trimmed_Ancestral/SYM+G4/fiveUTR_muscle_ires1_nt_hmm_sym_ancestral.fasta")
plotTree(ires1_trimmed_sym_tree, ftype = "i", fize = 0.5, lwd = 1)
labelnodes(1:ires1_trimmed_sym_tree$Nnode+Ntip(ires1_trimmed_sym_tree),
           1:ires1_trimmed_sym_tree$Nnode+Ntip(ires1_trimmed_sym_tree),
           interactive = FALSE, cex=0.3, shape = "ellipse")

## IRES 2 GTR - Untrimmed ## 

ires2_gtr_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES2_MUSCLE_Ancestral/GTR+I+G4/5UTR_MUSCLE_ires2_NT.fasta.raxml.ancestralTree")
ires2_gtr_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES2_MUSCLE_Ancestral/GTR+I+G4/fiveUTR_muscle_ires2_nt_gtr_ancestral.fasta")
plotTree(ires2_gtr_tree, ftype = "i", fize = 0.3, lwd = 1)
labelnodes(1:ires2_gtr_tree$Nnode+Ntip(ires2_gtr_tree),
           1:ires2_gtr_tree$Nnode+Ntip(ires2_gtr_tree),
           interactive = FALSE, cex=0.3, shape = "ellipse")

## IRES 2 GTR - Trimmed ## 

ires2_trimmed_gtr_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES2_MUSCLE_Trimmed_Ancestral/GTR+I+G4/5UTR_MUSCLE_ires2_NT_hmm.fasta.raxml.ancestralTree")
ires2_trimmed_gtr_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES2_MUSCLE_Trimmed_Ancestral/GTR+I+G4/fiveUTR_muscle_ires2_nt_hmm_gtr_ancestral.fasta")
plotTree(ires2_trimmed_gtr_tree, ftype = "i", fize = 0.5, lwd = 1)
labelnodes(1:ires2_trimmed_gtr_tree$Nnode+Ntip(ires2_trimmed_gtr_tree),
           1:ires2_trimmed_gtr_tree$Nnode+Ntip(ires2_trimmed_gtr_tree),
           interactive = FALSE, cex=0.3, shape = "ellipse")

## IRES 3 HKY - Trimmed ## 

ires3_trimmed_hky_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES3_MUSCLE_Trimmed_Ancestral/HKY+F+G4/5UTR_MUSCLE_ires3_NT_hmm.fasta.raxml.ancestralTree")
ires3_trimmed_hky_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES3_MUSCLE_Trimmed_Ancestral/HKY+F+G4/fiveUTR_muscle_ires3_nt_hmm_hky_ancestral.fasta")
plotTree(ires3_trimmed_hky_tree, ftype = "i", fize = 0.3, lwd = 1)
labelnodes(1:ires3_trimmed_hky_tree$Nnode+Ntip(ires3_trimmed_hky_tree),
           1:ires3_trimmed_hky_tree$Nnode+Ntip(ires3_trimmed_hky_tree),
           interactive = FALSE, cex=0.5, shape = "ellipse")

## IRES 3 GTR - Untrimmed ## 

ires3_gtr_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES3_MUSCLE_Ancestral/GTR+I+G4/5UTR_MUSCLE_ires3_NT.fasta.raxml.ancestralTree")
ires3_gtr_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES3_MUSCLE_Ancestral/GTR+I+G4/fiveUTR_muscle_ires3_nt_gtr_ancestral.fasta")
plotTree(ires3_gtr_tree, ftype = "i", fize = 0.3, lwd = 1)
labelnodes(1:ires3_gtr_tree$Nnode+Ntip(ires3_gtr_tree),
           1:ires3_gtr_tree$Nnode+Ntip(ires3_gtr_tree),
           interactive = FALSE, cex=0.5, shape = "ellipse")

## IRES 3 GTR - Trimmed ## 

ires3_trimmed_gtr_tree <- read.tree("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES3_MUSCLE_Trimmed_Ancestral/GTR+I+G4/5UTR_MUSCLE_ires3_NT_hmm.fasta.raxml.ancestralTree")
ires3_trimmed_gtr_seq <- read.fasta("~/Desktop/ViralAncestral/results/raxml-ng/data/5UTR_IRES3_MUSCLE_Trimmed_Ancestral/GTR+I+G4/fiveUTR_muscle_ires3_nt_hmm_gtr_ancestral.fasta")
plotTree(ires3_trimmed_gtr_tree, ftype = "i", fize = 0.3, lwd = 1)
labelnodes(1:ires3_trimmed_gtr_tree$Nnode+Ntip(ires3_trimmed_gtr_tree),
           1:ires3_trimmed_gtr_tree$Nnode+Ntip(ires3_trimmed_gtr_tree),
           interactive = FALSE, cex=0.5, shape = "ellipse")

