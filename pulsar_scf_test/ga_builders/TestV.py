import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_T=[-61.5805952694 ,-7.41082185633 ,-0.0144738837457 ,0.0 ,0.0 ,-1.23168557214 ,-1.23168557214 ,
-7.41082185633 ,-10.0090711421 ,-0.176890834734 ,0.0 ,0.0 ,-2.97722685358 ,-2.97722685358 ,
-0.0144738837457 ,-0.176890834734 ,-9.9440433417 ,0.0 ,0.0 ,-1.47179333871 ,-1.47179333871 ,
0.0 ,0.0 ,0.0 ,-9.87587599509 ,0.0 ,0.0 ,0.0 ,
0.0 ,0.0 ,0.0 ,0.0 ,-9.98754993509 ,-1.82223691348 ,1.82223691348 ,
-1.23168557214 ,-2.97722685358 ,-1.47179333871 ,0.0 ,-1.82223691348 ,-5.3002032523 ,-1.06717108047 ,
-1.23168557214 ,-2.97722685358 ,-1.47179333871 ,0.0 ,1.82223691348 ,-1.06717108047 ,-5.3002032523
]


def run(mm):
    tester = psr.PyTester("Testing Building of the nuclear-electron attraction matrix")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_GAV","V_INTS_KEY","LIBINT_V")
    bs=wf.system.get_basis_set("PRIMARY")
    T_maker=mm.get_module("PSR_GAV",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,wf,bs,bs)
        dims=_T[0].sizes()
        T=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                T[i,j]=_T[0].get_value([i,j])
        tester.test_double_vector("V",T.flatten(),corr_T)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
