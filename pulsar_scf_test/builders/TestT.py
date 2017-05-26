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
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    bs=wf.system.get_basis_set("PRIMARY")
    T_maker=mm.get_module("PSR_T",0)
    T=T_maker.calculate("???",0,wf,bs,bs)
    T=np.array(T[0].get_matrix())

    tester.test_double_vector("Kinetic energy of electrons",T.flatten(),corr_T)
    T_maker.options().change("FORCE_CACHE",True)
    T = T_maker.calculate("",0,wf,bs,bs)
    T=np.array(T[0].get_matrix())
    tester.test_double_vector("Kinetic energy cache works",T.flatten(),corr_T)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
