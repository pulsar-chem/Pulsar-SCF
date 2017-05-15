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

corr_J_DF=[18.0229,      3.57363, -1.72456e-17,    0.0108954,            0,     0.581893,     0.581893,
           3.57363,      8.77827, -3.32715e-17,    0.0282485,            0,      2.69973,      2.69973,
      -1.72456e-17, -3.32715e-17,      8.76241,  4.52858e-17,            0,      1.67385,     -1.67385,
         0.0108954,    0.0282485,  4.52858e-17,      8.77462, -1.33638e-51,      1.31232,      1.31232,
                 0,            0,            0, -1.33638e-51,      8.79376,  7.70372e-34,  7.70372e-34,
          0.581893,      2.69973,      1.67385,      1.31232,  7.70372e-34,      4.54872,     0.982587,
          0.581893,      2.69973,     -1.67385,      1.31232,  7.70372e-34,     0.982587,      4.54872
]

corr_K_DF=[-17.1411,     -3.23774,  1.24312e-17,    0.0054105,  -5.3926e-33,    -0.581427,    -0.581427,
           -3.23774,     -5.31048,   4.0495e-17,   -0.0366569, -3.38964e-32,       -1.486,       -1.486,
        2.18274e-17,  1.25361e-17,     -3.82877, -9.32084e-17,  9.93714e-18,     -0.59825,      0.59825,
          0.0054105,   -0.0366569,  -7.1537e-17,     -3.81804,  6.03972e-31,    -0.471044,    -0.471044,
       -9.62965e-35, -6.08594e-32,  9.93714e-18,  5.91838e-31,     -3.82841, -1.47209e-16,  1.47209e-16,
          -0.581427,       -1.486,     -0.59825,    -0.471044, -1.47209e-16,    -0.590544,    -0.385865,
          -0.581427,       -1.486,      0.59825,    -0.471044,  1.47209e-16,    -0.385865,    -0.590544
]

def run(mm):
    tester = psr.PyTester("Testing JK builder")
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_V","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    mm.change_option("PSR_JK","ERI_KEY","LIBINT_ERI")
    mm.change_option("PSR_Metric","METRIC_INTS_KEY","LIBINT_Metric")
    mm.change_option("PSR_3C2E","DF_INTS_KEY","LIBINT_3C2E")

    wf=make_wf()
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    JK=mm.get_module("PSR_JK",0).calculate("???",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())
    tester.test_double_vector("J",J.flatten(),corr_J)
    tester.test_double_vector("K",K.flatten(),corr_K)

    JK=mm.get_module("PSR_DFJK",0).calculate("",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())
    tester.test_double_vector("DF-J",J.flatten(),corr_J_DF)
    tester.test_double_vector("DF-K",K.flatten(),corr_K_DF)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
