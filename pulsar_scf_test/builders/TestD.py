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
    mm.load_module("pulsar_libint","NuclearElectron","V builder")
    mm.load_module("pulsar_libint","Kinetic","T builder")
    mm.load_module("pulsar_libint","Overlap","S builder")
    mm.load_module("pulsar_scf","TElectronic","T")
    mm.load_module("pulsar_scf","NuclearElectronic","V")
    mm.load_module("pulsar_scf","HCore","H")
    mm.load_module("pulsar_scf","Overlap","S")
    mm.load_module("pulsar_scf","CoreDensity","DCore")
    mm.change_option("V","V_INTS_KEY","V builder")
    mm.change_option("T","T_INTS_KEY","T builder")
    mm.change_option("S","S_INTS_KEY","S builder")
    mm.change_option("H","H_KEYS",["T","V"])
    mm.change_option("DCore","H_KEY","H")
    mm.change_option("DCore","S_KEY","S")
    wf=make_wf()
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("DCore",0).deriv(0,wf)[0]
    temp=guess.opdm.get(psr.Irrep.A,psr.Spin.alpha).get_matrix()
    D=np.array(temp)
    tester.test_double_vector("Core Guess",D.flatten(),corr_D)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
