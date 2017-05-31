import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_S=[ 1.0 ,0.236703936511 ,8.60874855256e-18 ,0.0 ,0.0 ,0.0384055905135 ,0.0384055905135 ,
0.236703936511 ,1.0 ,-2.46744544308e-17 ,0.0 ,0.0 ,0.386138781331 ,0.386138781331 ,
8.60874855256e-18 ,-2.46744544308e-17 ,1.0 ,0.0 ,0.0 ,0.209727649423 ,0.209727649423 ,
0.0 ,0.0 ,0.0 ,1.0 ,0.0 ,0.0 ,0.0 ,
0.0 ,0.0 ,0.0 ,0.0 ,1.0 ,0.268437641268 ,-0.268437641268 ,
0.0384055905135 ,0.386138781331 ,0.209727649423 ,0.0 ,0.268437641268 ,1.0 ,0.181760866822 ,
0.0384055905135 ,0.386138781331 ,0.209727649423 ,0.0 ,-0.268437641268 ,0.181760866822 ,1.0
]


def run(mm):
    tester = psr.PyTester("Testing Building of the overlap matrix")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_GAS","S_INTS_KEY","LIBINT_S")
    bs=wf.system.get_basis_set("PRIMARY")
    T_maker=mm.get_module("PSR_GAS",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,wf,bs,bs)
        dims=_T[0].sizes()
        T=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                T[i,j]=_T[0].get_value([i,j])
        tester.test_double_vector("S",T.flatten(),corr_S)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
