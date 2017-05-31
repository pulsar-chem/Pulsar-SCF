import pulsar as psr
import numpy as np
import os,sys
import pulsar_scf
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_J=[18.6894946288 ,3.56671424403 ,0.00363692677339 ,2.74532501597e-19 ,7.8062556419e-18 ,0.581617276796 ,0.581617276796 ,
3.56671424403 ,8.84604117875 ,-0.00195023073374 ,2.7668424914e-18 ,7.2858385991e-17 ,2.70860563168 ,2.70860563168 ,
0.00363692677339 ,-0.00195023073374 ,8.83741266219 ,-7.09056295686e-18 ,0.0 ,1.3111279185 ,1.3111279185 ,
2.74532501597e-19 ,2.7668424914e-18 ,-7.09056295686e-18 ,8.84116623159 ,-2.06941807663e-17 ,-2.02010179732e-18 ,2.22114814761e-18 ,
7.8062556419e-18 ,7.2858385991e-17 ,0.0 ,-2.06941807663e-17 ,8.83615337881 ,1.68060851305 ,-1.68060851305 ,
0.581617276796 ,2.70860563168 ,1.3111279185 ,-2.02010179732e-18 ,1.68060851305 ,4.55654075651 ,0.981027039215 ,
0.581617276796 ,2.70860563168 ,1.3111279185 ,2.22114814761e-18 ,-1.68060851305 ,0.981027039215 ,4.55654075651
]

corr_K=[-19.7014754062 ,-3.44227599747 ,-0.00276827584297 ,-2.40635548643e-19 ,-5.77337656849e-18 ,-0.592989598624 ,-0.592989598624 ,
-3.44227599747 ,-5.74400327672 ,-0.00811490273085 ,-2.45257499339e-17 ,-6.42931888284e-17 ,-1.5238076784 ,-1.5238076784 ,
-0.00276827584297 ,-0.00811490273085 ,-4.73190521254 ,1.18956106481e-16 ,-1.6826817717e-16 ,-0.548496139545 ,-0.548496139545 ,
-2.40635548643e-19 ,-2.45257499339e-17 ,1.18956106481e-16 ,-4.73965740936 ,3.05603779843e-16 ,-1.49883935721e-16 ,1.15249623178e-16 ,
-5.77337656849e-18 ,-6.42931888284e-17 ,-1.6826817717e-16 ,3.05603779843e-16 ,-4.73347899993 ,-0.705401791195 ,0.705401791195 ,
-0.592989598624 ,-1.5238076784 ,-0.548496139545 ,-1.49883935721e-16 ,-0.705401791195 ,-0.645615123776 ,-0.378177571713 ,
-0.592989598624 ,-1.5238076784 ,-0.548496139545 ,1.15249623178e-16 ,0.705401791195 ,-0.378177571713 ,-0.645615123776
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
    mm.change_option("PSR_GA_JK","ERI_KEY","LIBINT_ERI")
    mm.change_option("PSR_Sieve","ERI_INTS_KEY","LIBINT_ERI")
    mm.change_option("PSR_Metric","METRIC_INTS_KEY","LIBINT_Metric")
    mm.change_option("PSR_3C2E","DF_INTS_KEY","LIBINT_3C2E")

    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    JKmod=mm.get_module("PSR_GA_JK",0)
    for x in [False,True]:
        JKmod.options().change("FORCE_CACHE",x)
        JK=JKmod.calculate("???",0,guess,bs,bs)
        dims=JK[0].sizes()
        J=np.ndarray(dims)
        K=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                J[i,j]=JK[0].get_value([i,j])
                K[i,j]=JK[1].get_value([i,j])
        tester.test_double_vector("J",J.flatten(),corr_J)
        tester.test_double_vector("K",K.flatten(),corr_K)

    #JKmod=mm.get_module("PSR_DFJK",0)
    #JK=JKmod.calculate("",0,guess,bs,bs)
    #J=np.array(JK[0].get_matrix())
    #K=np.array(JK[1].get_matrix())
    #tester.test_double_vector("DF-J",J.flatten(),corr_J_DF)
    #tester.test_double_vector("DF-K",K.flatten(),corr_K_DF)

    #JKmod.options().change("FORCE_CACHE",True)
    #JK=JKmod.calculate("",0,guess,bs,bs)
    #J=np.array(JK[0].get_matrix())
    #K=np.array(JK[1].get_matrix())
    #tester.test_double_vector("Cache DF-J",J.flatten(),corr_J_DF)
    #tester.test_double_vector("Cache DF-K",K.flatten(),corr_K_DF)
    #tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
