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
    mm.change_option("PSR_Metric","METRIC_INTS_KEY","LIBINT_Metric")
    mm.change_option("PSR_3C2E","DF_INTS_KEY","LIBINT_3C2E")
    mm.change_option("PSR_G","JK_KEY","PSR_DFJK")
    wf=make_wf()
    wf.system=psr.apply_single_basis("PRIMARY","aug-cc-pvdz",wf.system)
    wf.system=psr.apply_single_basis("FITTING","aug-cc-pvdz-jkfit",wf.system)
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    #TODO test wf
    egy=mm.get_module("PSR_SCF",0).deriv(0,guess)[1][0]

    tester.test_double("Energy",-76.00333902407489,egy)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
