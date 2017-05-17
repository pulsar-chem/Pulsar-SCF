import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_X=[-0.346435643646409,                     0,     0.732191701095393,                     0,    -0.616437569630711,                     0,     0.176690380604319,
1.04111383320037,  1.04905412472594e-17,    -0.294748839341688,  2.43396912664122e-32,    -0.309011164496319, -1.67672878034322e-17,     0.444031298505881,
-4.76300102392361e-17,     0.859903586771966, -5.20096473450534e-17, -5.10280044091563e-16,  -4.4654750904489e-17,     -0.68873177878647, -1.19423209808818e-17,
0.48674950878724, -7.77227101784309e-17,     0.663796885540692,  1.47911418285228e-31,     0.623849682067528, -1.24139017149093e-16,        0.211348311338,
              -0, -2.15023909599824e-15,                     0,    -0.999999990239352,                     0, -5.37062348066963e-16,                     0,
-0.656546904528143,    -0.770785767203304,     -0.18081083589841,  1.13797858913338e-15,     0.148884012084651,    -0.384181281039476,     0.383787483246592,
-0.656546904528143,     0.770785767203303,     -0.18081083589841, -1.13797858913338e-15,     0.148884012084651,     0.384181281039476 ,    0.383787483246592,
]
corr_Xinv=[-0.150430044678067,                   0, 0.648535207242568,                    0, -0.678145783634254,                     0,     0.311273512061523,
            0.452074731095322, 5.4421242963472e-18, -0.26107233846135, 2.43396917415545e-32, -0.339944592320399, -2.1788659118268e-17,      0.78224508475469,
-2.06820074657735e-17,     0.446087774861987, -4.60672899857719e-17, -5.10280054052891e-16, -4.91248952643467e-17,    -0.894989227108336, -2.10386563273622e-17,
0.211357438810845, -4.03198118639011e-17,    0.587954834897832,  1.47911421172651e-31,     0.686299882353223, -1.61315458984653e-16,     0.372330009780075,
               -0, -1.11546851125005e-15,                     0,     -1.00000000976065,                     0, -6.97898703980758e-16,                     0,
-0.285087236237792,    -0.399856580524059,    -0.160152310871112,  1.13797861134819e-15,     0.163787981167732,     -0.49923369064347,     0.676114213953418,
-0.285087236237792,     0.399856580524059,    -0.160152310871112, -1.13797861134819e-15,     0.163787981167732,      0.49923369064347,     0.676114213953418
]


def run(mm):
    tester = psr.PyTester("Testing Building of the orthogonalizing matrix")
    wf=make_wf()
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_S","S_INTS_KEY","LIBINT_S")
    bs=wf.system.get_basis_set("PRIMARY")
    X_maker=mm.get_module("PSR_X",0)
    X=X_maker.calculate("???",0,wf,bs,bs)
    X,Xinv=np.array(X[0].get_matrix()),np.array(X[1].get_matrix())
    #Only defined up to a sign...
    for i,j in zip(X.flatten(),corr_X):
        desc="X: "+str(i)+","+str(j)
        tester.test_double(desc,abs(float(i)),abs(j))
    for i,j in zip(Xinv.flatten(),corr_Xinv):
        desc="Xinv: "+str(i)+","+str(j)
        tester.test_double(desc,abs(float(i)),abs(j))
    X_maker.options().change("FORCE_CACHE",True)
    X=X_maker.calculate("???",0,wf,bs,bs)
    X,Xinv=np.array(X[0].get_matrix()),np.array(X[1].get_matrix())
    #Only defined up to a sign...
    for i,j in zip(X.flatten(),corr_X):
        desc="Cache X: "+str(i)+","+str(j)
        tester.test_double(desc,abs(float(i)),abs(j))
    for i,j in zip(Xinv.flatten(),corr_Xinv):
        desc="Cache Xinv: "+str(i)+","+str(j)
        tester.test_double(desc,abs(float(i)),abs(j))
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
