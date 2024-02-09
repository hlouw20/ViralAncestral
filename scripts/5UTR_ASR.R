## Reconstruct Various Ancestral Sequences using RAxML 

cd: /Users/haileylouw/Desktop/ViralAncestral/data/raxml-ng

### IRES 1 

Trimmed:
- MFP tree: SYM+G4

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires1_NT_hmm.fasta --tree ./data/5UTR_MUSCLE_ires1_MFP_NT_hmm.fasta.treefile --msa-format FASTA --data-type DNA --model SYM+G4

- GTR + I + G4 

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires1_NT_hmm.fasta --tree ./data/5UTR_MUSCLE_ires1_GTR_NT_hmm.fasta.treefile --msa-format FASTA --data-type DNA --model GTR+I+G4

Untrimmed: 
- MFP tree: SYM+G4

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires1_NT.fasta --tree ./data/5UTR_MUSCLE_ires1_MFP_NT.fasta.treefile --msa-format FASTA --data-type DNA --model SYM+G4

- GTR + I + G4 

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires1_NT.fasta --tree ./data/5UTR_MUSCLE_ires1_GTR_NT.fasta.treefile --msa-format FASTA --data-type DNA --model GTR+I+G4

### IRES 2

Trimmed:
- MFP tree: TIM3e+I+G4
# not supported by RAxML
./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires2_NT_hmm.fasta --tree ./data/5UTR_MUSCLE_ires2_MFP_NT_hmm.fasta.treefile --msa-format FASTA --data-type DNA --model TIM3e+I+G4

- GTR + I + G4 

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires2_NT_hmm.fasta --tree ./data/5UTR_MUSCLE_ires2_GTR_NT_hmm.fasta.treefile --msa-format FASTA --data-type DNA --model GTR+I+G4


Untrimmed: 
- MFP tree: TPM3u+F+R4
# not supported by RAxML
./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires2_NT.fasta --tree ./data/5UTR_MUSCLE_ires2_MFP_NT.fasta.treefile --msa-format FASTA --data-type DNA --model TPM3u+F+R4

- GTR + I + G4 

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires2_NT.fasta --tree ./data/5UTR_MUSCLE_ires2_GTR_NT.fasta.treefile --msa-format FASTA --data-type DNA --model GTR+I+G4

### IRES 3 

Trimmed:
- MFP tree: HKY+F+G4

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires3_NT_hmm.fasta --tree ./data/5UTR_MUSCLE_ires3_MFP_NT_hmm.fasta.treefile --msa-format FASTA --data-type DNA --model HKY+F+G4

- GTR + I + G4 

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires3_NT_hmm.fasta --tree ./data/5UTR_MUSCLE_ires3_GTR_NT_hmm.fasta.treefile --msa-format FASTA --data-type DNA --model GTR+I+G4

Untrimmed: 
- MFP tree: TN+F+G4
# not supported by RAxML
./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires3_NT.fasta --tree ./data/5UTR_MUSCLE_ires3_MFP_NT.fasta.treefile --msa-format FASTA --data-type DNA --model TN+F+G4

- GTR + I + G4 

./raxml-ng --ancestral --msa ./data/5UTR_MUSCLE_ires3_NT.fasta --tree ./data/5UTR_MUSCLE_ires3_GTR_NT.fasta.treefile --msa-format FASTA --data-type DNA --model GTR+I+G4


