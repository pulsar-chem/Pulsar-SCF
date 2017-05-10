import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_T=[    29.0032,  -0.168011,       0,5.61892e-17,      0,-0.00841639,-0.00841639,
          -0.168011,   0.808128,       0,-2.2648e-17,      0,  0.0705173,  0.0705173,
                  0,          0, 2.52873,          0,      0,   0.147091,  -0.147091,
        5.61892e-17,-2.2648e-17,        0,   2.52873,      0,    0.11492,    0.11492,
                  0,          0,        0,         0,2.52873,          0,          0,
        -0.00841639,  0.0705173, 0.147091,   0.11492,      0,   0.760032,-0.00397974,
        -0.00841639,  0.0705173,-0.147091,   0.11492,      0,-0.00397974,   0.760032]


def run(mm):
    tester = psr.PyTester("Testing Building of the kinetic energy matrix")
    wf=make_wf()
    mm.load_module("pulsar_libint","Kinetic","T builder")
    mm.load_module("pulsar_scf","TElectronic","T")
    mm.change_option("T","T_INTS_KEY","T builder")
    bs=wf.system.get_basis_set("PRIMARY")
    T=mm.get_module("T",0).calculate("???",0,wf,bs,bs)
    T=np.array(T[0].get_matrix())
    tester.test_double_vector("Kinetic energy of electrons",T.flatten(),corr_T)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
