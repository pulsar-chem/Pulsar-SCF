import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_J=[18.6895,      3.56671,  1.73472e-17,   0.00363693,  3.85186e-33,     0.581617,     0.581617,
        3.56671,      8.84604,   1.4138e-16,  -0.00195023,   1.2326e-32,      2.70861,      2.70861,
    1.73472e-17,   1.4138e-16,      8.83615,  1.73472e-17, -2.92922e-17,      1.68061,     -1.68061,
     0.00363693,  -0.00195023,  1.73472e-17,      8.83741,   7.2415e-32,      1.31113,      1.31113,
    3.85186e-33,   1.2326e-32, -2.92922e-17,   7.2415e-32,      8.84117,  -4.7956e-18,   4.7956e-18,
       0.581617,      2.70861,      1.68061,      1.31113,  -4.7956e-18,      4.55654,     0.981027,
       0.581617,      2.70861,     -1.68061,      1.31113,   4.7956e-18,     0.981027,      4.55654
]

corr_K=[-19.7015,      -3.44228,  -2.23346e-17,   -0.00276828,   -2.6963e-33,      -0.59299,      -0.59299,
        -3.44228,        -5.744,  -5.35162e-16,    -0.0081149,  -9.24446e-33,      -1.52381,      -1.52381,
    -2.23346e-17,  -5.35162e-16,      -4.73348,  -5.44269e-16,  5.4817e-16,     -0.705402,    0.705402,
     -0.00276828,    -0.0081149,  -5.44269e-16,      -4.73191,  -1.24608e-30,     -0.548496,     -0.548496,
     -2.6963e-33,  -9.24446e-33,  5.4817e-16,  -1.24608e-30,      -4.73966, 2.21447e-16,  -2.21447e-16,
        -0.59299,      -1.52381,     -0.705402,     -0.548496, 2.21447e-16,     -0.645615,     -0.378178,
        -0.59299,      -1.52381,    0.705402,     -0.548496,  -2.21447e-16,     -0.378178,     -0.645615
]

def run(mm):
    tester = psr.PyTester("Testing JK builder")
    mm.load_module("pulsar_libint","NuclearElectron","V builder")
    mm.load_module("pulsar_libint","Kinetic","T builder")
    mm.load_module("pulsar_libint","Overlap","S builder")
    mm.load_module("pulsar_libint","ERI","ERI")
    mm.load_module("pulsar_libint","Metric","Metric builder")
    mm.load_module("pulsar_libint","DF3C2E","DF Ints builder")
    mm.load_module("pulsar_scf","TElectronic","T")
    mm.load_module("pulsar_scf","NuclearElectronic","V")
    mm.load_module("pulsar_scf","HCore","H")
    mm.load_module("pulsar_scf","Overlap","S")
    mm.load_module("pulsar_scf","CoreDensity","DCore")
    mm.load_module("pulsar_scf","JK","JK")
    mm.load_module("pulsar_scf","DFJK","DFJK")
    mm.load_module("pulsar_scf","Metric","Metric")
    mm.load_module("pulsar_scf","DFInts","DFInts")
    mm.change_option("V","V_INTS_KEY","V builder")
    mm.change_option("T","T_INTS_KEY","T builder")
    mm.change_option("S","S_INTS_KEY","S builder")
    mm.change_option("H","H_KEYS",["T","V"])
    mm.change_option("DCore","H_KEY","H")
    mm.change_option("DCore","S_KEY","S")
    mm.change_option("JK","ERI_KEY","ERI")
    mm.change_option("Metric","METRIC_INTS_KEY","Metric builder")
    mm.change_option("DFInts","DF_INTS_KEY","DF Ints builder")
    mm.change_option("DFJK","METRIC_KEY","Metric")
    mm.change_option("DFJK","DF_INTS_KEY","DFInts")
    wf=make_wf()
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("DCore",0).deriv(0,wf)[0]
    JK=mm.get_module("JK",0).calculate("???",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())
    tester.test_double_vector("J",J.flatten(),corr_J)
    tester.test_double_vector("K",K.flatten(),corr_K)

    JK=mm.get_module("DFJK",0).calculate("",0,guess,bs,bs)

    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
