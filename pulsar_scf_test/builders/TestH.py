import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_H=[-32.5774,   -7.57883,          0, -0.0144739,          0,    -1.2401,    -1.2401,
        -7.57883,   -9.20094,          0,  -0.176891,          0,   -2.90671,   -2.90671,
               0,          0,   -7.45882,          0,          0,   -1.67515,    1.67515,
      -0.0144739,  -0.176891,          0,   -7.41531,          0,   -1.35687,   -1.35687,
               0,          0,          0,          0,   -7.34714,          0,          0,
         -1.2401,   -2.90671,   -1.67515,   -1.35687,          0,   -4.54017,   -1.07115,
         -1.2401,   -2.90671,    1.67515,   -1.35687,          0,   -1.07115,   -4.54017
]


def run(mm):
    tester = psr.PyTester("Testing Building of the nuclear-electron attraction energy matrix")
    mm.load_module("pulsar_libint","NuclearElectron","V builder")
    mm.load_module("pulsar_scf","NuclearElectronic","V")
    mm.load_module("pulsar_libint","Kinetic","T builder")
    mm.load_module("pulsar_scf","TElectronic","T")
    mm.load_module("pulsar_scf","HCore","H")
    mm.change_option("V","V_INTS_KEY","V builder")
    mm.change_option("T","T_INTS_KEY","T builder")
    mm.change_option("H","H_KEYS",["T","V"])
    wf=make_wf()
    bs=wf.system.get_basis_set("PRIMARY")
    H=mm.get_module("H",0).calculate("???",0,wf,bs,bs)
    H=np.array(H[0].get_matrix())
    tester.test_double_vector("core Hamiltonian",H.flatten(),corr_H)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
