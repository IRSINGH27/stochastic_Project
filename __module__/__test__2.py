import numpy as np
def dnorm(x,mu,sd)->np.float64:
    a=1/(sd*(np.sqrt(2*np.pi)))
    b=-1*((x-mu)**2)/(2*(sd**2))
    return a*np.exp(b)
class Simulation():
    def __init__(self,p=0.2,q=0.3) -> None:
        assert 1-(2*p)-q>=0, 'Sum of proability is not 1'
        self.subsitute_dict={
        'a':{'a':1-(2*p)-q,'t':p,'g':q,'c':p},
        't':{'a':p,'t':1-(2*p)-q,'g':p,'c':q},
        'g':{'a':q,'t':p,'g':1-(2*p)-q,'c':p},
        'c':{'a':p,'t':q,'g':p,'c':1-(2*p)-q}
        }
        # self.indel={'insert':{'a':((1/3)/2)/4,'t':((1/3)/2)/4,'g':((1/3)/2)/4,'c':((1/3)/2)/4},'del':{'a':((1/3)/2)/4,'t':((1/3)/2)/4,'g':((1/3)/2)/4,'c':((1/3)/2)/4,},'None':2/3}
        self.np_round=np.vectorize(lambda x:round(x))
        self._nt_selector=np.vectorize(lambda x,y:x[y])
        return None

    def worker(self,seq_input:str,each_gen:int,offspring:int,gen_time:int,percent_mutation:float,percent_population=5) -> np.array:
        population=np.full(shape=(each_gen,1),fill_value=seq_input,dtype=object)
        for time_point in range(gen_time):
            _temp=population[:,time_point]
            np.random.seed()
            choosen_one=np.random.choice(_temp,size=round((percent_population/100)*each_gen),replace=True)
            gc_=np.vectorize(lambda x:(x.count('g')+x.count('c'))/len(x))(choosen_one)
            gc_=np.mean(gc_)
            childs=np.array([],dtype=object)
            for each_one in choosen_one:
                seq=self.reproduce(offspring=offspring,each_one=each_one,percent_mutation=percent_mutation)
                childs=np.append(childs,seq)
            p_=self.fitness_function(childrens=childs,gc_parent=gc_)
            np.random.seed()
            childs=np.random.choice(childs,size=each_gen,replace=True,p=p_)
            childs=np.reshape(childs,(each_gen,1))
            population=np.append(population,childs,axis=1)
        return population

    def reproduce(self,offspring,each_one,percent_mutation) ->np.array:
        results=np.array([],dtype=object)
        for _ in range(offspring):
            seq=each_one[3:-3]
            pos_number=round((percent_mutation/100)*len(seq))
            np.random.seed()
            pos=self.np_round(np.random.uniform(size=pos_number,low=0,high=len(seq)-1))
            _nt=self._nt_selector(seq,pos)
            seq=np.array([list(seq)])[0]
            # mutation
            for i in range(len(_nt)):
                sum_=0
                _proability=np.array([],dtype=np.float32)
                _mt=np.array([],dtype=str)
                for j in self.subsitute_dict[_nt[i]].keys():
                    if j!=_nt[i]:
                        sum_=sum_+self.subsitute_dict[_nt[i]][j]
                        _proability=np.append(_proability,self.subsitute_dict[_nt[i]][j])
                        _mt=np.append(_mt,j)
                    else:
                        sum_2=self.subsitute_dict[_nt[i]][j]
                np.random.seed()
                np.random.seed()
                t1=(1/(sum_))*np.log((1/np.random.rand()))
                # t11=(1/(self.indel['insert']['a']*8))*np.log(1/np.random.rand())
                t2=(1/sum_2)*np.log((1/np.random.rand()))
                t_happen=min([t1,t2])
                if t_happen==t1:
                    _proability=_proability/_proability.sum()
                    seq[pos[i]]=np.random.choice(_mt,p=_proability)
                # if t_happen==t11:
                #     event=np.random.choice(['insert','del'])
                #     if event=='insert':
                #         seq=np.insert(seq,pos[i],np.random.choice(tuple(self.indel[event].keys())))
                #         pos+=1
                #     elif event=='del':
                #         seq=np.delete(seq,pos[i])
                #         pos-=1
                #     else:
                #         continue    
                else:
                    continue
            seq=''.join(seq)
            seq=each_one[0:3]+seq+each_one[-3:]
            results=np.append(results,seq)
        return results
                    
    def fitness_function(self,childrens,gc_parent) ->np.array:
        diff_=np.vectorize(lambda x: ((x.count('g')+x.count('c'))/len(x))-gc_parent)(childrens)
        sd=np.std(diff_)
        if sd!=0:
            p_=np.vectorize(dnorm)(diff_,0,sd)
        else:
            p_=np.vectorize(dnorm)(diff_,0,1)
        p_=np.vectorize(lambda x:x/sum(p_))(p_)
        return p_

def main(p:float,q:float,i:int):
    import pandas as pd
    from Bio.Seq import Seq
    seq='atgcctaagtacctgccccctgacgccctcgtcgctctcatcaacaaggagttcggggccaacacgctcgtgcgcgcgaaggatgctgtcggcctcgtgaagccgcgcctgtctacaggttcctttgctctcgaccttcagctcggcggtggcttccccgaaggtgccatcactctgctcgaaggcgacaagggctcgtcaaagagctggaccatgaacaccatggccgcgatgttcctccagacgcacaagaacggtgtgttcatcctggtgaatgccgaaggcaccaacgaccacctgttcctcgaatcgctcggcgtcgataccgcgcgcaccttcttcctccagcccgagtcaggcgagcaggcctgggacgctgccatcaaagctgcgcagttcgctgagaaggtcttcatcggcgtcgattcgctcgatgcctgtgtgccgctcacggaacttgaaggagacgtgggcgatgccaagtacgcccctgccgccaagatgaacaacaagggcttccgcaagctcatctcggccatgaagcctgacctgaccagcacggatcagcgcgtcactgccgtgttcatcacccagctccgcgaagccatcggcgtcatgttcggtgatccgaagcgcagcgtcggtggcatgggcaaggcgttcgccgccatgaccatcatccgcctgtcgcgcatcaaggtgctgcgcaccgagggtgacaccgtcgctgaaaagaagagctacggcctggagatcgaggcgcacatcaccaagaacaagggatggggcgaaggcgaaaaggtgaagtggaccctctacaaagagaatcatgagggcttccgccgtggccagatcgacaacgtcaccgagctgattccgttcctgctcgtctacaagatcgcagacaagaagggtgcgtggatcaccctcggcaccgaccagtaccagggcgacaaggacctcgccgcccagctccgcatcaacgatgagctgcgggcgtggtgcatcgcccaggtgaaggaggcccacgccaagcgctacgagatgcaggaggaagtccctgccccgacgccgtccatcgtcaacaaaggcacctcggcgctgaagcgcctgcccaagaaaggcaagtaa'
    x=Simulation(p=p,q=q)
    y=x.worker(seq_input=seq,each_gen=100,offspring=1000,gen_time=100,percent_mutation=5,percent_population=5)
    results=pd.DataFrame.from_records(y)
    results.columns=[f'Gen_{i}' for i in results.columns]
    results.index=[f'org_{i}' for i in results.index]
    results_protein=results.applymap(lambda x:Seq(x).translate())
    results_protein=results_protein.applymap(lambda x:x[:x.find('*')])
    results2=results_protein.applymap(lambda x:len(x))
    results_protein=results_protein.astype(str)
    results=results.astype(str)
    results.to_parquet(f'{i}.result1.parquet')
    results_protein.to_parquet(f'{i}.result_protein1.parquet')
    results2.to_parquet(f'{i}_result21.parquet')