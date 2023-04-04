from __test__2 import main
def multi_function():
    from multiprocessing import Process
    _p=[]
    for p in range(10):
      _p.append(Process(target=main,args=(0.01,0.02,p)))
    for i in _p:
        i.start()
    for i in _p:
        i.join()
    return None

if __name__=='__main__':
    x=multi_function()