## Code for aligning 5UTR' Region Nucleotides through MACSE with default params
Directory: `cd ./Desktop/ViralAncestral/data`

- IRES type 1: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ./5UTR_Region/IRES_Subtypes/5UTR_ires1_NT.fasta -out_NT /Users/haileylouw/Desktop/ViralAncestral/data/5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires1_NT.fasta -out_AA /Users/haileylouw/Desktop/ViralAncestral/data/5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires1_AA.fasta`

- IRES type 2: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ./5UTR_Region/IRES_Subtypes/5UTR_ires2_NT.fasta -out_NT /Users/haileylouw/Desktop/ViralAncestral/data/5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires2_NT.fasta -out_AA /Users/haileylouw/Desktop/ViralAncestral/data/5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires2_AA.fasta`

- IRES type 3: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ./5UTR_Region/IRES_Subtypes/5UTR_ires3_NT.fasta -out_NT /Users/haileylouw/Desktop/ViralAncestral/data/5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires3_NT.fasta -out_AA /Users/haileylouw/Desktop/ViralAncestral/data/5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires3_AA.fasta`

## Code for using trimAl for aligned 5UTR' Region with gap threshold = 1 (the minimum fraction of sequences without a gap that you require to consider a column of “enough quality”)
Directory: `cd /Users/haileylouw/Desktop/ViralAncestral/data`

- IRES type 1 NT: `./Trimmed_Sequences/trimAl/source/trimal -in ./5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires1_NT.fasta -out ./Trimmed_Sequences/trimAl/5UTR_Region/FASTA_files/5UTR_MACSE_default_trimAL_ires1_NT.fasta -htmlout ./Trimmed_Sequences/trimAl/5UTR_Region/HTML_files/5UTR_MACSE_default_trimAL_ires1_NT.html -gt 1`

- IRES type 1 AA: `./Trimmed_Sequences/trimAl/source/trimal -in ./5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires1_AA.fasta -out ./Trimmed_Sequences/trimAl/5UTR_Region/FASTA_files/5UTR_MACSE_default_trimAL_ires1_AA.fasta -htmlout ./Trimmed_Sequences/trimAl/5UTR_Region/HTML_files/5UTR_MACSE_default_trimAL_ires1_AA.html -gt 1`

- IRES type 2 NT: `./Trimmed_Sequences/trimAl/source/trimal -in ./5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires2_NT.fasta -out ./Trimmed_Sequences/trimAl/5UTR_Region/FASTA_files/5UTR_MACSE_default_trimAL_ires2_NT.fasta -htmlout ./Trimmed_Sequences/trimAl/5UTR_Region/HTML_files/5UTR_MACSE_default_trimAL_ires2_NT.html -gt 1`

- IRES type 2 AA: `./Trimmed_Sequences/trimAl/source/trimal -in ./5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires2_AA.fasta -out ./Trimmed_Sequences/trimAl/5UTR_Region/FASTA_files/5UTR_MACSE_default_trimAL_ires2_AA.fasta -htmlout ./Trimmed_Sequences/trimAl/5UTR_Region/HTML_files/5UTR_MACSE_default_trimAL_ires2_AA.html -gt 1`

- IRES type 3 NT: `./Trimmed_Sequences/trimAl/source/trimal -in ./5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires3_NT.fasta -out ./Trimmed_Sequences/trimAl/5UTR_Region/FASTA_files/5UTR_MACSE_default_trimAL_ires3_NT.fasta -htmlout ./Trimmed_Sequences/trimAl/5UTR_Region/HTML_files/5UTR_MACSE_default_trimAL_ires3_NT.html -gt 1`

- IRES type 3 AA: `./Trimmed_Sequences/trimAl/source/trimal -in ./5UTR_Region/Aligned_Sequences/MACSE_Default/5UTR_MACSE_default_ires3_AA.fasta -out ./Trimmed_Sequences/trimAl/5UTR_Region/FASTA_files/5UTR_MACSE_default_trimAL_ires3_AA.fasta -htmlout ./Trimmed_Sequences/trimAl/5UTR_Region/HTML_files/5UTR_MACSE_default_trimAL_ires3_AA.html -gt 1`