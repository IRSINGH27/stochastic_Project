!#/bin/bash
# musc='/home/irsingh/miniconda3/envs/scanners/bin/muscle'
# for i in {0..100}
# do
# $musc -super5 test_fasta_2_codon.fasta -output test_fasta_2_codon.@.mfa -perm all  -threads 4 -perturb $i
# done
iqtree='/home/irsingh/miniconda3/envs/env_phylo/bin/iqtree'
file='/mnt/c/Users/inder/Documents/GitHub/stochastic_Project/__module__/test_fasta_2_codon.msa'
$iqtree -s $file --seqtype DNA --safe --runs 1 -T 4 -B 20000 -m MFP --ancestral --prefix test_tree_2_codon -redo
