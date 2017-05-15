import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_D=[1.06501,    -0.285217, -2.19377e-18,   -0.0195534, -9.13187e-33,    0.0334495,    0.0334495,
      -0.285217,      1.24897, -1.14244e-16,      0.11356,  7.43117e-32,    -0.144281,    -0.144281,
   -2.19377e-18, -1.14244e-16,      1.12587, -3.24462e-16, -1.24947e-16,    -0.146131,     0.146131,
     -0.0195534,      0.11356, -3.24462e-16,      1.06606,  4.06834e-31,   -0.0993586,   -0.0993586,
   -9.13187e-33,  7.43117e-32, -1.24947e-16,  4.06834e-31,            1, -1.00919e-16,  1.00919e-16,
      0.0334495,    -0.144281,    -0.146131,   -0.0993586, -1.00919e-16,    0.0426801,   0.00474611,
      0.0334495,    -0.144281,     0.146131,   -0.0993586,  1.00919e-16,   0.00474611,    0.0426801
]


def run(mm):
    tester = psr.PyTester("Testing core guess")
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_V","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")

    wf=make_wf()
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
