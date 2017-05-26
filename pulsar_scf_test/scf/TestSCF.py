import pulsar as psr
import numpy as np
import os,sys
import functools
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf

def run(mm):
    tester = psr.PyTester("Testing pulsar SCF")
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_V","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    mm.change_option("PSR_JK","ERI_KEY","LIBINT_ERI")
    mm.change_option("PSR_Sieve","ERI_INTS_KEY","LIBINT_ERI")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    wf.system=psr.apply_single_basis("PRIMARY","aug-cc-pvdz",wf.system)
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    #TODO test wf
    egy=mm.get_module("PSR_SCF",0).deriv(0,guess)[1][0]
    tester.test_double("RHF Energy",-76.00335405845436,egy)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
