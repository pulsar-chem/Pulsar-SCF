import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
corr_metric=[1.04643709 ,3.42919963053 ,6.32625272634e-18 ,0.0 ,0.0 ,2.60526240572 ,2.60526240572 ,
3.42919963053 ,26.4352252164 ,-1.8004034429e-16 ,0.0 ,0.0 ,25.3420821293 ,25.3420821293 ,
6.32625272634e-18 ,-1.8004034429e-16 ,5.78479783655 ,0.0 ,0.0 ,3.2924421174 ,3.2924421174 ,
0.0 ,0.0 ,0.0 ,5.78479783655 ,0.0 ,0.0 ,0.0 ,
0.0 ,0.0 ,0.0 ,0.0 ,5.78479783655 ,4.21411005387 ,-4.21411005387 ,
2.60526240572 ,25.3420821293 ,3.2924421174 ,0.0 ,4.21411005387 ,39.9325707859 ,26.6712894369 ,
2.60526240572 ,25.3420821293 ,3.2924421174 ,0.0 ,-4.21411005387 ,26.6712894369 ,39.9325707859
]


def run(mm):
    tester = psr.PyTester("Testing Building of the density fitting metric matrix")
    wf=psr.make_wf("sto-3g","""
    O 0.0 -0.07579 0.0
    H 0.86681 0.60144 0.0
    H -0.86681 0.60144 0.0
    """)
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_Metric","METRIC_INTS_KEY","LIBINT_Metric")
    bs=wf.system.get_basis_set("PRIMARY")
    metric_builder = mm.get_module("PSR_Metric",0)
    metric=metric_builder.calculate("???",0,wf,bs,bs)
    metric=np.array(metric[0].get_matrix())
    tester.test_double_vector("Metric matrix",metric.flatten(),corr_metric)
    metric_builder.options().change("FORCE_CACHE",True)
    metric=metric_builder.calculate("???",0,wf,bs,bs)
    metric=np.array(metric[0].get_matrix())
    tester.test_double_vector("Metric matrix cache works",metric.flatten(),corr_metric)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
