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

## Revision Fall 2024 

Goal: reconstruct Triticum Mosaic Virus 

Method: reconstruct both replicase gene and 5'UTR (region)

  - 4 Final Ancestral Sequences:
  
    1. 5'UTR Nucleotides using GTR + I + G4 model of evolution
      a. Results in: 
       
    2. 5'UTR Nucleotides using SYM + G4 model of evolution 
      a. Results in: 
      
    3. Replicse Gene AA using model finder in IQ-Tree 
      a. Results in: 
      
    4. Replicase Gene NT using model finder in IQ-Tree
    
Analysis Steps: 

  - Aquire original sequences for both 5'UTR and Replicase in Google Drive from Helena and Aurelie 
  
  - Split the 2 respective areas of the genome among the 5 IRES types 
    a. By splitting up the files for these regions by IRES type, we were not left with much data (sequences) in each subset by IRES type.  
      i. IRES types 1 - 3 were the most well-studied / understood at the time of this analysis, so we decided to focus on mainly these 3 IRES types
      ii. IRES type 1 had 22 sequences in its subset, while IRES types 2 and 3 had fewer than 10 sequences in their subsets.
      iii. Because we were not confident in the downstream analyses which would be conducted using very few sequences, we decided to focus on only IRES type 1.  
      
  - Attach mitovirus sequences (2 sequences) as outgroup to the IRES subset
    
  - Align using MUSCLE bioinformatics software 
    a. We would like to trim these sequences using HmmCleaner, but do to computational restrictions in the dependencies required for this 
    
  - Build Tree using IQ-Tree
    a. For NT, we find that GTR+I+G4 is widely used in the literature and considered the most biologically reasonable model of evolution across many organism types 
      i. A limitation to this is that viral evolution is one of the least well-understood evolutionary patterns in biology, so it is unclear whether this model's feasibility is applicable to viruses as well. 
      
    b. For NT, we also allow IQ-Tree's ModelFinderPlus software choose the best-fitting model for the data, from a computational perspective.  
      i. For the 5-UTR MSA, this was SYM+G4 
      
    c. For the Replicase Gene NT, we constructed 2 trees using GTR+I+F+G4 and ModelFinder, but ModelFinder found that GTR+I+F+G4 was the best-fitting model as well.  
    
    c. All trees were constructed using 100 bootstrap replicates (not ultrafast)
      
  - Root tree on Mitovirus outgroup     
      
  - Reconstruct Ancestral Sequence for Triticum using RAxML-NG
    a. We perform four (4) reconstructions, supplying each respective MSA and corresponding tree (built under a specific evolutionary pattern) to the software and running the program.
      i. Note that, for some evolutionary models used in IQ-Tree, there are not yet suppoerted by RAxML-NG.  Though RAxML-NG is more widely-used for ASR and has shown to have the most feasibly downstream analyses, IQ-Tree can also perform ASR, and if a viable option if the model of evolution cannot be used in RAxML-NG.  
      
  - Validate trees + sequences with biologists
  
  - Send sequences to be synthesized in a laboratory 
    a. If the sequence is viable (doesn't die off immediately), we consider this a "success"
    
    
## To-do
- ASR for Replicase NT tree
- ASR for Replicase AA tree(s)
- Send sequences and trees to Helena for verification   
  
  
  
  
  




