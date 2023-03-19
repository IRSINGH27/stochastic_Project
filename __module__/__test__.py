
from Bio.Seq import Seq
import numpy as np
seq_id='''atgcctaagtacctgccccctgacgccctcgtcgctctcatcaacaaggagttcggggcc
aacacgctcgtgcgcgcgaaggatgctgtcggcctcgtgaagccgcgcctgtctacaggt
tcctttgctctcgaccttcagctcggcggtggcttccccgaaggtgccatcactctgctc
gaaggcgacaagggctcgtcaaagagctggaccatgaacaccatggccgcgatgttcctc
cagacgcacaagaacggtgtgttcatcctggtgaatgccgaaggcaccaacgaccacctg
ttcctcgaatcgctcggcgtcgataccgcgcgcaccttcttcctccagcccgagtcaggc
gagcaggcctgggacgctgccatcaaagctgcgcagttcgctgagaaggtcttcatcggc
gtcgattcgctcgatgcctgtgtgccgctcacggaacttgaaggagacgtgggcgatgcc
aagtacgcccctgccgccaagatgaacaacaagggcttccgcaagctcatctcggccatg
aagcctgacctgaccagcacggatcagcgcgtcactgccgtgttcatcacccagctccgc
gaagccatcggcgtcatgttcggtgatccgaagcgcagcgtcggtggcatgggcaaggcg
ttcgccgccatgaccatcatccgcctgtcgcgcatcaaggtgctgcgcaccgagggtgac
accgtcgctgaaaagaagagctacggcctggagatcgaggcgcacatcaccaagaacaag
ggatggggcgaaggcgaaaaggtgaagtggaccctctacaaagagaatcatgagggcttc
cgccgtggccagatcgacaacgtcaccgagctgattccgttcctgctcgtctacaagatc
gcagacaagaagggtgcgtggatcaccctcggcaccgaccagtaccagggcgacaaggac
ctcgccgcccagctccgcatcaacgatgagctgcgggcgtggtgcatcgcccaggtgaag
gaggcccacgccaagcgctacgagatgcaggaggaagtccctgccccgacgccgtccatc
gtcaacaaaggcacctcggcgctgaagcgcctgcccaagaaaggcaagtaa'''.replace('\n','')
seq_0=Seq(seq_id,)
seq_0.translate()
p=0.2
p_dict={'a':{'a':1-3*p,'t':p,'g':p,'c':p},
        't':{'a':p,'t':1-3*p,'g':p,'c':p},
        'g':{'a':p,'t':p,'g':1-3*p,'c':p},
        'c':{'a':p,'t':p,'g':p,'c':1-3*p}}

_temp_seq=seq_0[3:-3]
pos=round(np.random.uniform(low=0,high=len(_temp_seq)))
_nt=_temp_seq[pos]
change=np.random.choice(p_dict[_nt.lower()],p)