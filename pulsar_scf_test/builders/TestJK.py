import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
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

corr_J_DF=[18.0229053915 ,3.57363169646 ,0.0108954065104 ,4.82060117179e-19 ,-6.24500451352e-17 ,0.581892530675 ,0.581892530675 ,
3.57363169646 ,8.77826703769 ,0.0282485396006 ,3.09973915447e-18 ,-1.80411241502e-16 ,2.69972828925 ,2.69972828925 ,
0.0108954065104 ,0.0282485396006 ,8.7746211517 ,-2.66069793265e-36 ,1.04083408559e-16 ,1.31231557633 ,1.31231557633 ,
4.82060117179e-19 ,3.09973915447e-18 ,-2.66069793265e-36 ,8.79376127459 ,0.0 ,8.55677617915e-19 ,8.55677617915e-19 ,
-6.24500451352e-17 ,-1.80411241502e-16 ,1.04083408559e-16 ,0.0 ,8.76240530959 ,1.67384832211 ,-1.67384832211 ,
0.581892530675 ,2.69972828925 ,1.31231557633 ,8.55677617915e-19 ,1.67384832211 ,4.54872424809 ,0.982587263459 ,
0.581892530675 ,2.69972828925 ,1.31231557633 ,8.55677617915e-19 ,-1.67384832211 ,0.982587263459 ,4.54872424809
]

corr_K_DF=[-17.1411078652 ,-3.23774192706 ,0.0054104957526 ,-8.06657548717e-21 ,6.59194920871e-17 ,-0.581427105593 ,-0.581427105593 ,
-3.23774192706 ,-5.31047761713 ,-0.0366569284149 ,-2.4276068061e-17 ,1.73472347598e-16 ,-1.48600490082 ,-1.48600490082 ,
0.0054104957526 ,-0.0366569284149 ,-3.81804037541 ,1.13252485978e-16 ,-2.25514051877e-16 ,-0.471044224452 ,-0.471044224452 ,
-8.06657548717e-21 ,-2.4276068061e-17 ,1.13252485978e-16 ,-3.82840822358 ,2.83722449278e-16 ,-1.53195795396e-16 ,1.17538154177e-16 ,
6.59194920871e-17 ,1.73472347598e-16 ,-2.25514051877e-16 ,2.83722449278e-16 ,-3.82876949891 ,-0.598250198707 ,0.598250198707 ,
-0.581427105593 ,-1.48600490082 ,-0.471044224452 ,-1.53195795396e-16 ,-0.598250198707 ,-0.590544173381 ,-0.385864677735 ,
-0.581427105593 ,-1.48600490082 ,-0.471044224452 ,1.17538154177e-16 ,0.598250198707 ,-0.385864677735 ,-0.590544173381
]

def run(mm):
    tester = psr.PyTester("Testing JK builder")
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_V","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_T","T_INTS_KEY","LIBINT_T")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    mm.change_option("PSR_JK","ERI_KEY","LIBINT_ERI")
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
    JKmod=mm.get_module("PSR_JK",0)
    JK=JKmod.calculate("???",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())

    tester.test_double_vector("J",J.flatten(),corr_J)
    tester.test_double_vector("K",K.flatten(),corr_K)

    JKmod.options().change("FORCE_CACHE",True)
    JK=JKmod.calculate("???",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())
    tester.test_double_vector("Cache J",J.flatten(),corr_J)
    tester.test_double_vector("Cache K",K.flatten(),corr_K)

    #I'm lazy and am going to use the STO-3G for the fitting basis too
    guess.system=psr.apply_single_basis("FITTING","sto-3g",guess.system)

    JKmod=mm.get_module("PSR_DFJK",0)
    JK=JKmod.calculate("",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())

    tester.test_double_vector("DF-J",J.flatten(),corr_J_DF)
    tester.test_double_vector("DF-K",K.flatten(),corr_K_DF)

    JKmod.options().change("FORCE_CACHE",True)
    JK=JKmod.calculate("",0,guess,bs,bs)
    J=np.array(JK[0].get_matrix())
    K=np.array(JK[1].get_matrix())
    tester.test_double_vector("Cache DF-J",J.flatten(),corr_J_DF)
    tester.test_double_vector("Cache DF-K",K.flatten(),corr_K_DF)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
