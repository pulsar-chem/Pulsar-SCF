import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

corr_F=[
-18.8132695464909752, -4.8726875509324366, -0.0115290259221462, 0.0000000000000000, 0.0000000000000000, -0.8067320791902532, -0.8067320791902530,
-4.8726875509324366, -1.7909028272382308, -0.1808697909925568, -0.0000000000000000, 0.0000000000000001, -0.5790558023646937, -0.5790558023646937,
-0.0115290259221462, -0.1808697909925568, 0.2391242158958038, 0.0000000000000000, -0.0000000000000001, -0.1828690746748323, -0.1828690746748314,
0.0000000000000000, -0.0000000000000000, 0.0000000000000000, 0.3091070826153004, -0.0000000000000001, -0.0000000000000000, 0.0000000000000000,
0.0000000000000000, 0.0000000000000001, -0.0000000000000001, -0.0000000000000001, 0.1939648923276049, -0.1708882956899287, 0.1708882956899287,
-0.8067320791902532, -0.5790558023646937, -0.1828690746748323, -0.0000000000000000, -0.1708882956899287, -0.1450343928538382, -0.1846681705831821,
-0.8067320791902530, -0.5790558023646937, -0.1828690746748314, 0.0000000000000000, 0.1708882956899287, -0.1846681705831821, -0.1450343928538347
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
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    mm.change_option("PSR_GA_JK","ERI_KEY","LIBINT_ERI")
    mm.change_option("PSR_Sieve","ERI_INTS_KEY","LIBINT_ERI")
    mm.change_option("PSR_DCore","H_KEY","PSR_GAH")
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    T_maker=mm.get_module("PSR_GAF",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,guess,bs,bs)
        dims=_T[0].sizes()
        T=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                T[i,j]=_T[0].get_value([i,j])
        tester.test_double_vector("F",T.flatten(),corr_F)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
