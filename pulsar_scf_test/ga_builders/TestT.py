import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_T=[
29.0031999455 ,-0.168010939316 ,5.61892214316e-17 ,0.0 ,0.0 ,-0.00841638518545 ,-0.00841638518545 ,
-0.168010939316 ,0.80812795493 ,-2.26480264056e-17 ,0.0 ,0.0 ,0.070517338519 ,0.070517338519 ,
5.61892214316e-17 ,-2.26480264056e-17 ,2.52873119819 ,0.0 ,0.0 ,0.114920380257 ,0.114920380257 ,
0.0 ,0.0 ,0.0 ,2.52873119819 ,0.0 ,0.0 ,0.0 ,
0.0 ,0.0 ,0.0 ,0.0 ,2.52873119819 ,0.147090552413 ,-0.147090552413 ,
-0.00841638518545 ,0.070517338519 ,0.114920380257 ,0.0 ,0.147090552413 ,0.760031883567 ,-0.00397973672704 ,
-0.00841638518545 ,0.070517338519 ,0.114920380257 ,0.0 ,-0.147090552413 ,-0.00397973672704 ,0.760031883567 ,
]

def run(mm):
    tester = psr.PyTester("Testing Building of the kinetic energy matrix")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_GAT","T_INTS_KEY","LIBINT_T")
    bs=wf.system.get_basis_set("PRIMARY")
    T_maker=mm.get_module("PSR_GAT",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,wf,bs,bs)
        dims=_T[0].sizes()
        T=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                T[i,j]=_T[0].get_value([i,j])
        tester.test_double_vector("T",T.flatten(),corr_T)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
