import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_D=[1.06501171016 ,-0.285216559127 ,-0.0195534221141 ,-1.14044644195e-19 ,2.81104002801e-18 ,0.0334495398727 ,0.0334495398727 ,
-0.285216559127 ,1.24896555611 ,0.113559685833 ,-6.27961222406e-18 ,-2.63427531904e-17 ,-0.144280859156 ,-0.144280859156 ,
-0.0195534221141 ,0.113559685833 ,1.06606416639 ,-4.2105380444e-17 ,-2.96395838375e-16 ,-0.0993585909193 ,-0.0993585909193 ,
-1.14044644195e-19 ,-6.27961222406e-18 ,-4.2105380444e-17 ,1.0 ,-1.58785880025e-16 ,1.91563909624e-16 ,-1.49891497978e-16 ,
2.81104002801e-18 ,-2.63427531904e-17 ,-2.96395838375e-16 ,-1.58785880025e-16 ,1.12586952884 ,-0.146131347598 ,0.146131347598 ,
0.0334495398727 ,-0.144280859156 ,-0.0993585909193 ,1.91563909624e-16 ,-0.146131347598 ,0.0426801126591 ,0.00474610662419 ,
0.0334495398727 ,-0.144280859156 ,-0.0993585909193 ,-1.49891497978e-16 ,0.146131347598 ,0.00474610662419 ,0.0426801126591
]


def run(mm):
    tester = psr.PyTester("Testing core guess")
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_V","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")

    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    temp=guess.opdm.get(psr.Irrep.A,psr.Spin.alpha).get_matrix()
    D=np.array(temp)
    tester.test_double_vector("Core Guess",D.flatten(),corr_D)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
