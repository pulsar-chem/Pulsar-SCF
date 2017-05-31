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
    tester = psr.PyTester("Testing Building of the core Hamiltonian")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_GAV","V_INTS_KEY","LIBINT_V")
    mm.change_option("PSR_GAT","T_INTS_KEY","LIBINT_T")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    mm.change_option("PSR_GA_JK","ERI_KEY","LIBINT_ERI")
    mm.change_option("PSR_Sieve","ERI_INTS_KEY","LIBINT_ERI")
    mm.change_option("PSR_DCore","H_KEY","PSR_GAH")
    bs=wf.system.get_basis_set("PRIMARY")
    guess=mm.get_module("PSR_DCore",0).deriv(0,wf)[0]
    T_maker=mm.get_module("PSR_GAG",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,guess,bs,bs)
        dims=_T[0].sizes()
        T=np.ndarray(dims)
        for i in range(dims[0]):
            for j in range(dims[1]):
                T[i,j]=_T[0].get_value([i,j])
        tester.test_double_vector("G",T.flatten(),corr_G)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
