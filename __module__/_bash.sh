!#/bin/bash
# musc='/home/irsingh/miniconda3/envs/scanners/bin/muscle'
# for i in {0..100}
# do
# $musc -super5 test_fasta_2.fasta -output test_fasta_2.@.mfa -perm all  -threads 4 -perturb $i
# done
iqtree='/home/irsingh/miniconda3/envs/env_phylo/bin/iqtree'
file='/mnt/c/Users/inder/Documents/GitHub/stochastic_Project/__module__/test_fasta_2.msa'
$iqtree -s $file --seqtype AA --safe --runs 1 -T 4 -B 5000 -m MFP+LMSS --ancestral --prefix test_tree_2
