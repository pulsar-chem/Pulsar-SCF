import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_T=[-61.5806,   -7.41082,          0, -0.0144739,          0,   -1.23169,   -1.23169,
        -7.41082,   -10.0091,          0,  -0.176891,          0,   -2.97723,   -2.97723,
               0,          0,   -9.98755,          0,          0,   -1.82224,    1.82224,
      -0.0144739,  -0.176891,          0,   -9.94404,          0,   -1.47179,   -1.47179,
               0,          0,          0,          0,   -9.87588,          0,          0,
        -1.23169,   -2.97723,   -1.82224,   -1.47179,          0,    -5.3002,   -1.06717,
        -1.23169,   -2.97723,    1.82224,   -1.47179,          0,   -1.06717,    -5.3002,
]


def run(mm):
    tester = psr.PyTester("Testing Building of the nuclear-electron attraction energy matrix")
    wf=make_wf()
    mm.load_module("pulsar_libint","NuclearElectron","V builder")
    mm.load_module("pulsar_scf","NuclearElectronic","V")
    mm.change_option("V","V_INTS_KEY","V builder")
    bs=wf.system.get_basis_set("PRIMARY")
    T=mm.get_module("V",0).calculate("???",0,wf,bs,bs)
    T=np.array(T[0].get_matrix())
    tester.test_double_vector("nuclear-electron attraction",T.flatten(),corr_T)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
