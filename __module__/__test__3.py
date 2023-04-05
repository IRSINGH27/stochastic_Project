import numpy as np
from Bio.Seq import Seq
def dnorm(x,mu,sd):
    a=1/(sd*(np.sqrt(2*np.pi)))
    b=-1*((x-mu)**2)/(2*(sd**2))
    return a*np.exp(b)
def len_seq(x):
    x=Seq(x).translate()
    x=x[0:x.find('*')]
    return len(x)
class Simulation():
    def __init__(self) -> None:
        # self.subsitute_dict={'a':{'a':1-(2*p)-q,'t':p,'g':q,'c':p},
        # 't':{'a':p,'t':1-(2*p)-q,'g':p,'c':q},
        # 'g':{'a':q,'t':p,'g':1-(2*p)-q,'c':p},
        # 'c':{'a':p,'t':q,'g':p,'c':1-(2*p)-q}}
        self.indel={'insert':{'a':((1/3)/2)/4,'t':((1/3)/2)/4,'g':((1/3)/2)/4,'c':((1/3)/2)/4},'del':{'a':((1/3)/2)/4,'t':((1/3)/2)/4,'g':((1/3)/2)/4,'c':((1/3)/2)/4,},'None':2/3}
        self.np_round=np.vectorize(lambda x:round(x))
        self._nt_selector=np.vectorize(lambda x,y:x[y])

    def worker(self,seq_input:str,each_gen:int,offspring:int,gen_time:int,percent_mutation:float,percent_population=5) -> np.array:
        population=np.full(shape=(each_gen,1),fill_value=seq_input,dtype=object)
        for time_point in range(gen_time):
            _temp=population[:,time_point]
            np.random.seed()
            choosen_one=np.random.choice(_temp,size=round((percent_population/100)*each_gen),replace=True)
            _len=np.vectorize(len_seq)(choosen_one)
            _len=np.mean(_len)
            childs=np.array([],dtype=object)
            for each_one in choosen_one:
                seq=self.reproduce(offspring=offspring,each_one=each_one,percent_mutation=percent_mutation)
                childs=np.append(childs,seq)
            p_=self.fitness_function(childrens=childs,len_parents=_len)
            np.random.seed()
            childs=np.random.choice(childs,size=each_gen,replace=True,p=p_)
            # childs=np.random.choice(childs,size=each_gen,replace=True)
            childs=np.reshape(childs,(each_gen,1))
            population=np.append(population,childs,axis=1)
        return population

    def reproduce(self,offspring,each_one,percent_mutation):
        results=np.array([],dtype=object)
        for _ in range(offspring):
            seq=each_one[3:-3]
            pos_number=round((percent_mutation/100)*len(seq))
            np.random.seed()
            pos=self.np_round(np.random.uniform(size=pos_number,low=0,high=len(seq)-1))
            seq=np.array([list(seq)])[0]
            # mutation
            for i in range(len(pos)):
                np.random.seed()
                t1=1/(((1/3)/2)/4) * np.log(1/np.random.uniform(0,1))
                np.random.seed()
                t2=1/(((1/3)/2)/4) * np.log(1/np.random.uniform(0,1))
                np.random.seed()
                t3=1/(1/3) * np.log(1/np.random.uniform(0,1))                
                t_happen=min([t1,t2,t3])
                if t_happen==t1:
                    seq=np.insert(seq,i,np.random.choice(['a','t','g','c']))
                elif t_happen==t2:
                    seq=np.delete(seq,i)
                else:
                    continue
            seq=''.join(seq)
            seq=each_one[0:3]+seq+each_one[-3:]
            results=np.append(results,seq)
        return results
                    
    def fitness_function(self,childrens,len_parents):
        childs_len=np.vectorize(len_seq)(childrens)
        diff_=childs_len-len_parents
        sd=np.std(diff_)
        if sd!=0:
            p_=np.vectorize(dnorm)(diff_,0,sd)
        else:
            p_=np.vectorize(dnorm)(diff_,0,1)
        p_=np.vectorize(lambda x:x/sum(p_))(p_)
        return p_
def main(i:int):
    import pandas as pd
    from Bio.Seq import Seq
    seq='atgcctaagtacctgccccctgacgccctcgtcgctctcatcaacaaggagttcggggccaacacgctcgtgcgcgcgaaggatgctgtcggcctcgtgaagccgcgcctgtctacaggttcctttgctctcgaccttcagctcggcggtggcttccccgaaggtgccatcactctgctcgaaggcgacaagggctcgtcaaagagctggaccatgaacaccatggccgcgatgttcctccagacgcacaagaacggtgtgttcatcctggtgaatgccgaaggcaccaacgaccacctgttcctcgaatcgctcggcgtcgataccgcgcgcaccttcttcctccagcccgagtcaggcgagcaggcctgggacgctgccatcaaagctgcgcagttcgctgagaaggtcttcatcggcgtcgattcgctcgatgcctgtgtgccgctcacggaacttgaaggagacgtgggcgatgccaagtacgcccctgccgccaagatgaacaacaagggcttccgcaagctcatctcggccatgaagcctgacctgaccagcacggatcagcgcgtcactgccgtgttcatcacccagctccgcgaagccatcggcgtcatgttcggtgatccgaagcgcagcgtcggtggcatgggcaaggcgttcgccgccatgaccatcatccgcctgtcgcgcatcaaggtgctgcgcaccgagggtgacaccgtcgctgaaaagaagagctacggcctggagatcgaggcgcacatcaccaagaacaagggatggggcgaaggcgaaaaggtgaagtggaccctctacaaagagaatcatgagggcttccgccgtggccagatcgacaacgtcaccgagctgattccgttcctgctcgtctacaagatcgcagacaagaagggtgcgtggatcaccctcggcaccgaccagtaccagggcgacaaggacctcgccgcccagctccgcatcaacgatgagctgcgggcgtggtgcatcgcccaggtgaaggaggcccacgccaagcgctacgagatgcaggaggaagtccctgccccgacgccgtccatcgtcaacaaaggcacctcggcgctgaagcgcctgcccaagaaaggcaagtaa'
    x=Simulation()
    y=x.worker(seq_input=seq,each_gen=10,offspring=100,gen_time=10,percent_mutation=1,percent_population=50)
    results=pd.DataFrame.from_records(y)
    results.columns=[f'Gen_{i}' for i in results.columns]
    results.index=[f'org_{i}' for i in results.index]
    results_protein=results.applymap(lambda x:Seq(x).translate())
    results_protein=results_protein.applymap(lambda x:x[:x.find('*')])
    results2=results_protein.applymap(lambda x:len(x))
    results_protein=results_protein.astype(str)
    results=results.astype(str)
    results.to_parquet(f'{i}.indel_result1.parquet')
    results_protein.to_parquet(f'{i}.indel_result_protein1.parquet')
    results2.to_parquet(f'{i}.indel_result21.parquet')
