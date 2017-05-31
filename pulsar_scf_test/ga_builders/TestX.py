import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_X=[-0.346435651004 ,0.0 ,0.732191695572 ,0.0 ,-0.616437587962 ,0.0 ,0.17669037509 ,
1.04111385331 ,2.6162235582e-33 ,-0.294748838856 ,-4.35212583698e-18 ,-0.309011170974 ,2.0954398593e-33 ,0.444031291464 ,
0.486749512201 ,2.38255910495e-32 ,0.663796906528 ,-1.55798642899e-16 ,0.623849674077 ,1.90828849469e-32 ,0.211348312475 ,
-0.0 ,-1.3501281658e-16 ,0.0 ,1.0 ,0.0 ,-1.08137256272e-16 ,0.0 ,
-0.0 ,0.859903605249 ,0.0 ,0.0 ,0.0 ,0.688731772918 ,0.0 ,
-0.656546895999 ,-0.770785758067 ,-0.180810835516 ,1.34116583553e-17 ,0.148884011101 ,0.384181290575 ,0.383787490884 ,
-0.656546895999 ,0.770785758067 ,-0.180810835516 ,1.34116583553e-17 ,0.148884011101 ,-0.384181290575 ,0.383787490884
]
corr_Xinv=[-0.150430046051 ,0.0 ,0.648535191315 ,0.0 ,-0.678145791828 ,0.0 ,0.311273500152 ,
0.452074734352 ,1.357204879e-33 ,-0.261072333589 ,-4.35212583698e-18 ,-0.339944593445 ,2.72297016324e-33 ,0.782245066832 ,
0.211357437734 ,1.23598796884e-32 ,0.587954843482 ,-1.55798642899e-16 ,0.686299861447 ,2.47977178196e-32 ,0.372330009157 ,
-0.0 ,-7.0039906496e-17 ,0.0 ,1.0 ,0.0 ,-1.40521581212e-16 ,0.0 ,
-0.0 ,0.446087783612 ,0.0 ,0.0 ,0.0 ,0.894989211841 ,0.0 ,
-0.285087229082 ,-0.399856575036 ,-0.160152307807 ,1.34116583553e-17 ,0.163787977194 ,0.499233698772 ,0.67611422264 ,
-0.285087229082 ,0.399856575036 ,-0.160152307807 ,1.34116583553e-17 ,0.163787977194 ,-0.499233698772 ,0.67611422264
]


def run(mm):
    tester = psr.PyTester("Testing Building of the orthogonalizer matrix")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_GAS","S_INTS_KEY","LIBINT_S")
    bs=wf.system.get_basis_set("PRIMARY")
    T_maker=mm.get_module("PSR_GAX",0)
    for x in [False,True]:
        T_maker.options().change("FORCE_CACHE",x)
        _T=T_maker.calculate("???",0,wf,bs,bs)
        dims=_T[0].sizes()
        for i in range(dims[0]):
            for j in range(dims[1]):
                x=_T[0].get_value([i,j])
                xcorr=corr_X[i*7+j]
                name="X "+str(x)+" "+str(xcorr)
                tester.test_double(name,abs(x),abs(xcorr))
                xinv=_T[1].get_value([i,j])
                xinvcorr=corr_Xinv[i*7+j]
                name="Xinv "+str(xinv)+" "+str(xinvcorr)
                tester.test_double(name,abs(xinv),abs(xinvcorr))
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
