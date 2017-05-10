import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_S=[2.18747927145885,   0.36996404054683,  0.205903188897061, 0.0606886002263033, 0.0606886002263033,
        0.36996404054683,  0.903994671141006,  0.559166380724417,  0.322205724679542,  0.322205724679542,
       0.205903188897061,  0.559166380724417,    1.5683964676283,  0.337647462105432,  0.337647462105432,
      0.0606886002263033,  0.322205724679542,  0.337647462105432,  0.880116997794134,  0.133647040945866,
      0.0606886002263033,  0.322205724679542,  0.337647462105432,  0.133647040945866,  0.880116997794134
]

def run(mm):
    tester = psr.PyTester("Testing Schwarz Screen")
    wf=make_wf()
    mm.load_module("pulsar_libint","ERI","ERI builder")
    mm.load_module("pulsar_scf","SchwarzScreen","Sieve")
    mm.change_option("Sieve","ERI_INTS_KEY","ERI builder")
    bs=wf.system.get_basis_set("PRIMARY")
    T=mm.get_module("Sieve",0).calculate("???",0,wf,bs,bs)
    S=np.array(T[0].get_matrix())
    tester.test_double_vector("Schwarz Screening",S.flatten(),corr_S)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
