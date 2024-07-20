# Viral Ancestral
Infer a phylogenetic tree from viral sequences, and use the tree to estimate the ancestral protein sequence of this group.

## Overall Workflow 

- Accessed sequences from original project's [Google Drive] (https://drive.google.com/drive/u/2/folders/0AHmGEazpc0PkUk9PVA) 

- Assessed both 5'UTR and Replicase sequences separately 
  1. 5'UTR AA: `data / picornaviridae_5UTR_with_outgroup.fasta`
  2. Replicase AA: `data / picornaviridae_replicase_aminoacid_with_outgroup.fasta`
  3. Replicase NT: `data / picornaviridae_replicase_nucleotide_with_outgroup.fasta`

- Split up dataset into IRES types 1-5

- Reattached mitovirus outgroup to IRES subsets (mitovirus is technically type 0 since it doesn't have an IRES)

  - Align sequences with outgroup using MUSCLE
  
  - Further split sequences such that they were either trimmed (HmmCleaner) or untrimmed 
  
    - Note that HmmCleaner's newest release didn't allow us to trim the sequences after reformatting the data, so our end product only includes untrimmed sequences 
    
    - In investigating the inferential advantages of trimming, it wasn't obvious that trimming would lead to more accurate downstream analyses, so we still hold a certain level of confidence in our untrimmed alignments 
    
- Construct trees using IQ-Tree with 2 renditions of each sequence / trimming: 
  (1) Run `-MFP` so that the best-fitting model is used to construct the tree 
  (2) Run `GTR+F+I+G4` additionally, after performing exhaustive literature searches as to what might be the most biologically-informed model to construct the tree with 
  
- Root trees on the mitovirus outgroup (using `phytools` in R)

- Use the rooted tree + sequence alignment in RAxML-NG to perform ASR 

  - This will result in 2 trees for each rendition of the alignment - one from `MFP` and one from `GTR`
  
- Obtain ancestral sequences, visualize, send to collaborators for confirmation in biological sense pertaining to the sequences 

- If sequences and tree have biological logic, send ancestral sequence(s) to be sequenced and tested for viability 






