import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_S=[  1,     0.236704,            0,  8.60875e-18,            0,    0.0384056,    0.0384056,
   0.236704,            1,            0, -2.46745e-17,            0,     0.386139,     0.386139,
          0,            0,            1,            0,            0,     0.268438,    -0.268438,
8.60875e-18, -2.46745e-17,            0,            1,            0,     0.209728,     0.209728,
          0,            0,            0,            0,            1,            0,            0,
  0.0384056,     0.386139,     0.268438,     0.209728,            0,            1,     0.181761,
  0.0384056,     0.386139,    -0.268438,     0.209728,            0,     0.181761,            1
]


def run(mm):
    tester = psr.PyTester("Testing Building of the overlap matrix")
    wf=make_wf()
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    bs=wf.system.get_basis_set("PRIMARY")
    S=mm.get_module("PSR_S",0).calculate("???",0,wf,bs,bs)
    S=np.array(S[0].get_matrix())
    tester.test_double_vector("Overlap matrix",S.flatten(),corr_S)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
