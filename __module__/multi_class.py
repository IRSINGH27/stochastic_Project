from __test__3 import main
def multi_function():
    from multiprocessing import Process
    _p=[]
    for p in range(10):
<<<<<<< HEAD
      _p.append(Process(target=main,args=(0.1,0.2,p)))
=======
      _p.append(Process(target=main,args=(p,)))
>>>>>>> a1235e8e088902c1c6f57413ec52f33bf6da3f1d
    for i in _p:
        i.start()
    for i in _p:
        i.join()
    return None

if __name__=='__main__':
    x=multi_function()