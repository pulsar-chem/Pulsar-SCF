import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

corr_G=[
13.7641257769471217,2.7061452446076224,0.0029448578099142,0.0000000000000000,0.0000000000000000,0.4333698765888307,0.4333698765888310,
2.7061452446076224,7.4100403594670343,-0.0039789564090952,-0.0000000000000000,0.0000000000000001,2.3276537099584802,2.3276537099584802,
0.0029448578099142,-0.0039789564090952,7.6544363589564721,0.0000000000000000,-0.0000000000000001,1.1740038827363963,1.1740038827363972,
0.0000000000000000,-0.0000000000000000,0.0000000000000000,7.6562518791485559,-0.0000000000000001,-0.0000000000000000,0.0000000000000000,
0.0000000000000000,0.0000000000000001,-0.0000000000000001,-0.0000000000000001,7.6527836287267510,1.5042580641319037,-1.5042580641319034,
0.4333698765888307,2.3276537099584802,1.1740038827363963,-0.0000000000000000,1.5042580641319037,4.3951369741154096,0.8864826449821217,
0.4333698765888310,2.3276537099584802,1.1740038827363972,0.0000000000000000,-1.5042580641319034,0.8864826449821217,4.3951369741154105
]

def run(mm):
    tester = psr.PyTester("Testing G builder")
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
    Fmod=mm.get_module("PSR_G",0)
    for x in [False,True]:
        Fmod.options().change("FORCE_CACHE",x)
        F=Fmod.calculate("???",0,guess,bs,bs)
        F=np.array(F[0].get_matrix())
        tester.test_double_vector("G",F.flatten(),corr_G)

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
