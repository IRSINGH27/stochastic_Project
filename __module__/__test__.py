
from Bio.Seq import Seq
import numpy as np
seq_id='''atgcctaagtacctgccccctgacgccctcgtcgctctcatcaacaaggagttcggggccaacacgctcgtgcgcgcgaaggatgctgtcggcctcgtgaagccgcgcctgtctacaggttcctttgctctcgaccttcagctcggcggtggcttccccgaaggtgccatcactctgctcgaaggcgacaagggctcgtcaaagagctggaccatgaacaccatggccgcgatgttcctccagacgcacaagaacggtgtgttcatcctggtgaatgccgaaggcaccaacgaccacctgttcctcgaatcgctcggcgtcgataccgcgcgcaccttcttcctccagcccgagtcaggcgagcaggcctgggacgctgccatcaaagctgcgcagttcgctgagaaggtcttcatcggcgtcgattcgctcgatgcctgtgtgccgctcacggaacttgaaggagacgtgggcgatgccaagtacgcccctgccgccaagatgaacaacaagggcttccgcaagctcatctcggccatgaagcctgacctgaccagcacggatcagcgcgtcactgccgtgttcatcacccagctccgcgaagccatcggcgtcatgttcggtgatccgaagcgcagcgtcggtggcatgggcaaggcgttcgccgccatgaccatcatccgcctgtcgcgcatcaaggtgctgcgcaccgagggtgacaccgtcgctgaaaagaagagctacggcctggagatcgaggcgcacatcaccaagaacaagggatggggcgaaggcgaaaaggtgaagtggaccctctacaaagagaatcatgagggcttccgccgtggccagatcgacaacgtcaccgagctgattccgttcctgctcgtctacaagatcgcagacaagaagggtgcgtggatcaccctcggcaccgaccagtaccagggcgacaaggacctcgccgcccagctccgcatcaacgatgagctgcgggcgtggtgcatcgcccaggtgaaggaggcccacgccaagcgctacgagatgcaggaggaagtccctgccccgacgccgtccatcgtcaacaaaggcacctcggcgctgaagcgcctgcccaagaaaggcaagtaa'''
seq_0=Seq(seq_id)
seq_0.translate()
p=0.2
q=0.3
subsitute_dict={'a':{'a':1-2*p-q,'t':p,'g':q,'c':p},
        't':{'a':p,'t':1-2*p-q,'g':p,'c':q},
        'g':{'a':q,'t':p,'g':1-2*p-q,'c':p},
        'c':{'a':p,'t':q,'g':p,'c':1-2*p-q}}

indel={'ins':{'a':0.25,'t':0.25,'g':0.25,'c':0.25},'del':{'a':0.25,'t':0.25,'g':0.25,'c':0.25}}


_temp_seq=seq_0[3:-3]
pos=round(np.random.uniform(low=0,high=len(_temp_seq)))
_nt=_temp_seq[pos]
change=np.random.choice(subsitute_dict[_nt.lower()],p)
from Bio import SeqIO
SeqIO.SeqRecord