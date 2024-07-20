## Code for aligning Replicase gene Nucleotides through MACSE with default params
Directory: `cd ./Desktop/ViralAncestral/data`

- IRES type 1: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ./Replicase_Gene/IRES_Subtypes/Replicase_ires1_NT.fasta -out_NT /Users/haileylouw/Desktop/ViralAncestral/data/Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires1_NT.fasta -out_AA /Users/haileylouw/Desktop/ViralAncestral/data/Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires1_AA.fasta`

- IRES type 2: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ./Replicase_Gene/IRES_Subtypes/Replicase_ires2_NT.fasta -out_NT /Users/haileylouw/Desktop/ViralAncestral/data/Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires2_NT.fasta -out_AA /Users/haileylouw/Desktop/ViralAncestral/data/Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires2_AA.fasta`

- IRES type 3: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ./Replicase_Gene/IRES_Subtypes/Replicase_ires3_NT.fasta -out_NT /Users/haileylouw/Desktop/ViralAncestral/data/Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires3_NT.fasta -out_AA /Users/haileylouw/Desktop/ViralAncestral/data/Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires3_AA.fasta`

## Code for using trimAl for aligned Replicase gene with gap threshold = 1 (the minimum fraction of sequences without a gap that you require to consider a column of “enough quality”)
Directory: `cd /Users/haileylouw/Desktop/ViralAncestral/data`

- IRES type 1 NT: `./Trimmed_Sequences/trimAl/source/trimal -in ./Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires1_NT.fasta -out ./Replicase_Gene/Trimmed_Sequences/trimAl/Replicase_Gene/FASTA_files/Rep_MACSE_default_trimAL_ires1_NT.fasta -htmlout ./Trimmed_Sequences/trimAl/Replicase_Gene/HTML_files/Rep_MACSE_default_trimAL_ires1_NT.html -gt 1`

- IRES type 1 AA: `./Trimmed_Sequences/trimAl/source/trimal -in ./Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires1_AA.fasta -out ./Trimmed_Sequences/trimAl/Replicase_Gene/FASTA_files/Rep_MACSE_default_trimAL_ires1_AA.fasta -htmlout ./Trimmed_Sequences/trimAl/Replicase_Gene/HTML_files/Rep_MACSE_default_trimAL_ires1_AA.html -gt 1`

- IRES type 2 NT: `./Trimmed_Sequences/trimAl/source/trimal -in ./Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires2_NT.fasta -out ./Replicase_Gene/Trimmed_Sequences/trimAl/Replicase_Gene/FASTA_files/Rep_MACSE_default_trimAL_ires2_NT.fasta -htmlout ./Trimmed_Sequences/trimAl/Replicase_Gene/HTML_files/Rep_MACSE_default_trimAL_ires2_NT.html -gt 1`

- IRES type 2 AA: `./Trimmed_Sequences/trimAl/source/trimal -in ./Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires2_AA.fasta -out ./Trimmed_Sequences/trimAl/Replicase_Gene/FASTA_files/Rep_MACSE_default_trimAL_ires2_AA.fasta -htmlout ./Trimmed_Sequences/trimAl/Replicase_Gene/HTML_files/Rep_MACSE_default_trimAL_ires2_AA.html -gt 1`

- IRES type 3 NT: `./Trimmed_Sequences/trimAl/source/trimal -in ./Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires3_NT.fasta -out ./Replicase_Gene/Trimmed_Sequences/trimAl/Replicase_Gene/FASTA_files/Rep_MACSE_default_trimAL_ires3_NT.fasta -htmlout ./Trimmed_Sequences/trimAl/Replicase_Gene/HTML_files/Rep_MACSE_default_trimAL_ires3_NT.html -gt 1`

- IRES type 3 AA: `./Trimmed_Sequences/trimAl/source/trimal -in ./Replicase_Gene/Aligned_Sequences/MACSE_Default/Rep_MACSE_default_ires3_AA.fasta -out ./Trimmed_Sequences/trimAl/Replicase_Gene/FASTA_files/Rep_MACSE_default_trimAL_ires3_AA.fasta -htmlout ./Trimmed_Sequences/trimAl/Replicase_Gene/HTML_files/Rep_MACSE_default_trimAL_ires3_AA.html -gt 1`