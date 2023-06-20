## Code for aligning Replicase gene through MACSE with default params
Directory: `cd/Users/haileylouw/Desktop/ViralAncestral/ReplicaseGene`

- IRES type 1: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ires1_for_repNT.fasta`
- IRES type 2: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ires2_for_repNT.fasta`
- IRES type 3: `java -jar ./macse_v2.07.jar -prog alignSequences -seq ires3_for_repNT.fasta`

## Code for using trimAl for aligned Replicase gene with gap threshold = 1 (the minimum fraction of sequences without a gap that you require to consider a column of “enough quality”)
Directory: `cd/Users/haileylouw/Desktop/ViralAncestral/ReplicaseGene/trimAL/source`

- IRES type 1 NT: `./trimal -in ires1_for_repNT_NT.fasta -out trimmedIres1NT.fasta -htmlout trimmedIres1NT.html -gt 1`

- IRES type 1 AA: `./trimal -in ires1_for_repNT_AA.fasta -out trimmedIres1AA.fasta -htmlout trimmedIres1AA.html -gt 1`

- IRES type 2 NT: `./trimal -in ires2_for_repNT_NT.fasta -out trimmedIres2NT.fasta -htmlout trimmedIres2NT.html -gt 1`

- IRES type 2 AA: `./trimal -in ires2_for_repNT_AA.fasta -out trimmedIres2AA.fasta -htmlout trimmedIres2AA.html -gt 1`

- IRES type 3 NT: `./trimal -in ires3_for_repNT_NT.fasta -out trimmedIres3NT.fasta -htmlout trimmedIres3NT.html -gt 1`

- IRES type 3 AA: `./trimal -in ires3_for_repNT_AA.fasta -out trimmedIres3AA.fasta -htmlout trimmedIres3AA.html -gt 1`