import pulsar as psr
import numpy as np
import os,sys
import functools
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf

def run(mm):
    tester = psr.PyTester("Testing pulsar SCF")
    mm.load_module("pulsar_libint","NuclearElectron","V builder")
    mm.load_module("pulsar_libint","Kinetic","T builder")
    mm.load_module("pulsar_libint","Overlap","S builder")
    mm.load_module("pulsar_libint","ERI","ERI")
    mm.load_module("pulsar_scf","TElectronic","T")
    mm.load_module("pulsar_scf","NuclearElectronic","V")
    mm.load_module("pulsar_scf","F","F")
    mm.load_module("pulsar_scf","G","G")
    mm.load_module("pulsar_scf","HCore","H")
    mm.load_module("pulsar_scf","Overlap","S")
    mm.load_module("pulsar_scf","Orthogonalizer","X")
    mm.load_module("pulsar_scf","CoreDensity","DCore")
    mm.load_module("pulsar_scf","JK","JK")
    mm.load_module("pulsar_scf","SCF","SCF")
    mm.change_option("V","V_INTS_KEY","V builder")
    mm.change_option("T","T_INTS_KEY","T builder")
    mm.change_option("S","S_INTS_KEY","S builder")
    mm.change_option("X","S_KEY","S")
    mm.change_option("H","H_KEYS",["T","V"])
    mm.change_option("G","JK_KEY","JK")
    mm.change_option("F","H_KEY","H")
    mm.change_option("F","G_KEY","G")
    mm.change_option("SCF","F_KEY","F")
    mm.change_option("SCF","S_KEY","S")
    mm.change_option("SCF","H_KEY","H")
    mm.change_option("SCF","X_KEY","X")
    mm.change_option("DCore","H_KEY","H")
    mm.change_option("DCore","S_KEY","S")
    mm.change_option("JK","ERI_KEY","ERI")
    wf=make_wf()
    wf.system=psr.apply_single_basis("PRIMARY","aug-cc-pvdz",wf.system)
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("DCore",0).deriv(0,wf)[0]
    egy=mm.get_module("SCF",0).deriv(0,guess)
    print(egy)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
