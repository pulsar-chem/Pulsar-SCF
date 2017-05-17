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
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_V","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    wf=make_wf()
    bs=wf.system.get_basis_set("PRIMARY")
    Hbuilder=mm.get_module("PSR_H",0)
    H=Hbuilder.calculate("???",0,wf,bs,bs)
    H=np.array(H[0].get_matrix())
    tester.test_double_vector("core Hamiltonian",H.flatten(),corr_H)
    Hbuilder.options().change("FORCE_CACHE",True)
    H=Hbuilder.calculate("???",0,wf,bs,bs)
    H=np.array(H[0].get_matrix())
    tester.test_double_vector("core Hamiltonian",H.flatten(),corr_H)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
