import pulsar as psr
import numpy as np
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_metric=[     1.04643710731328,      3.42919970498325,                     0,  6.32625284043065e-18,                     0,      2.60526241540254,      2.60526241540254,
3.42919970498325,      26.4352259268833,                     0, -1.80040348466822e-16,                     0,      25.3420823544302,      25.3420823544302,
               0,                     0,      5.78479794947718,                     0,                     0,      4.21411007580444,     -4.21411007580444,
6.32625284043065e-18, -1.80040348466822e-16,                     0,      5.78479794947718,                     0,      3.29244213453588,      3.29244213453588,
               0,                     0,                     0,                     0,      5.78479794947718,                     0,                     0,
2.60526241540254,      25.3420823544302,      4.21411007580444,      3.29244213453588,                     0,      39.9325704220625,      26.6712891938732,
2.60526241540254,      25.3420823544302,     -4.21411007580444,      3.29244213453588 ,                    0,     26.6712891938732,      39.9325704220625,
]


def run(mm):
    tester = psr.PyTester("Testing Building of the density fitting metric matrix")
    wf=make_wf()
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
