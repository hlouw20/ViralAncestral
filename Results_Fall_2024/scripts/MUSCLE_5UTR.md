## Code for aligning 5'UTR Region Nucleotides through MUSCLE 

Directory: `cd ./Desktop/ViralAncestral/Results_Fall_2024/data`

IRES type 1: `./muscle -in ./5UTR_Region/IRES_Subtypes/5UTR_ires1_NT.fasta -out ./5UTR_Region/Aligned_Sequences/MUSCLE/5UTR_MUSCLE_ires1_NT.fasta`

IRES type 2: `./muscle -in ./5UTR_Region/IRES_Subtypes/5UTR_ires2_NT.fasta -out ./5UTR_Region/Aligned_Sequences/MUSCLE/5UTR_MUSCLE_ires2_NT.fasta`

IRES type 3: `./muscle -in ./5UTR_Region/IRES_Subtypes/5UTR_ires3_NT.fasta -out ./5UTR_Region/Aligned_Sequences/MUSCLE/5UTR_MUSCLE_ires3_NT.fasta`

IRES type 1 with mitovirus outgroup: `./muscle -in ./5UTR_Region/IRES_Subtypes/IRES1_with_outgroup.fasta -out ./5UTR_Region/Aligned_Sequences/MUSCLE/IRES1_with_outgroup_aligned.fasta`
