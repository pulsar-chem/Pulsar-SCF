import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_H=[-32.5773953239 ,-7.57883279565 ,-0.0144738837457 ,0.0 ,0.0 ,-1.24010195733 ,-1.24010195733 ,
-7.57883279565 ,-9.20094318714 ,-0.176890834734 ,0.0 ,0.0 ,-2.90670951506 ,-2.90670951506 ,
-0.0144738837457 ,-0.176890834734 ,-7.4153121435 ,0.0 ,0.0 ,-1.35687295846 ,-1.35687295846 ,
0.0 ,0.0 ,0.0 ,-7.3471447969 ,0.0 ,0.0 ,0.0 ,
0.0 ,0.0 ,0.0 ,0.0 ,-7.45881873689 ,-1.67514636106 ,1.67514636106 ,
-1.24010195733 ,-2.90670951506 ,-1.35687295846 ,0.0 ,-1.67514636106 ,-4.54017136873 ,-1.0711508172 ,
-1.24010195733 ,-2.90670951506 ,-1.35687295846 ,0.0 ,1.67514636106 ,-1.0711508172 ,-4.54017136873
]

def run(mm):
    tester = psr.PyTester("Testing Building of the core Hamiltonian")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_GAV","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_GAT","T_INTS_KEY","LIBINT_T")
    bs=wf.system.get_basis_set("PRIMARY")
    T_maker=mm.get_module("PSR_GAH",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,wf,bs,bs)
        dims=_T[0].sizes()
        T=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                T[i,j]=_T[0].get_value([i,j])
        tester.test_double_vector("H",T.flatten(),corr_H)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
