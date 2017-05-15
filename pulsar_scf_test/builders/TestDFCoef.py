import pulsar as psr
import os,sys
sys.path.insert(0,os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pulsar_scf_test.TestCommon import make_wf
corr_ints=[2.0222080644043485,0.3698737616321260,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0603113001496220,0.0603113001496220,0.3698737616321261,0.6880699277230037,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.1953020171135259,0.1953020171135259,0.0000000000000000,0.0000000000000000,0.6880961673192155,0.0000000000000000,0.0000000000000000,0.1158094565282536,-0.1158094565282536,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.6880961673192155,0.0000000000000000,0.0904807723083827,0.0904807723083827,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.6880961673192155,0.0000000000000000,0.0000000000000000,0.0603113001496220,0.1953020171135259,0.1158094565282536,0.0904807723083827,0.0000000000000000,0.2981405659865650,0.0666249438292437,0.0603113001496220,0.1953020171135259,-0.1158094565282536,0.0904807723083827,0.0000000000000000,0.0666249438292437,0.2981405659865650,-0.3596205464209228,0.0026517383573205,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0001550172420762,0.0001550172420762,0.0026517383573204,0.5557326887499007,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.2205474476263449,0.2205474476263449,0.0000000000000000,0.0000000000000000,0.5553415758529688,0.0000000000000000,0.0000000000000000,0.1545524770562344,-0.1545524770562344,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.5553415758529688,0.0000000000000000,0.1207503074915997,0.1207503074915997,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.5553415758529688,0.0000000000000000,0.0000000000000000,0.0001550172420762,0.2205474476263448,0.1545524770562344,0.1207503074915997,0.0000000000000000,0.5065615886119672,0.0968814121700369,0.0001550172420762,0.2205474476263448,-0.1545524770562344,0.1207503074915997,0.0000000000000000,0.0968814121700369,0.5065615886119672,0.0000000000000000,0.0000000000000000,0.0599992976854208,0.0000000000000000,0.0000000000000000,0.0024777344020380,-0.0024777344020380,0.0000000000000000,0.0000000000000000,0.3858070096404518,0.0000000000000000,0.0000000000000000,0.0917620803284007,-0.0917620803284007,0.0599992976854208,0.3858070096404518,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.1647690501730896,0.1647690501730896,0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0455239678491316,-0.0455239678491316,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0024777344020380,0.0917620803284007,0.1647690501730896,0.0455239678491316,0.0000000000000000,0.3093857568235218,0.0000000000000000,-0.0024777344020380,-0.0917620803284007,0.1647690501730896,-0.0455239678491316,0.0000000000000000,-0.0000000000000000,-0.3093857568235218,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0599992976854208,0.0000000000000000,0.0019358291541309,0.0019358291541309,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.3858070096404518,0.0000000000000000,0.0716927973382896,0.0716927973382896,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0455239678491316,-0.0455239678491316,0.0599992976854208,0.3858070096404518,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.1420687775652751,0.1420687775652751,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0019358291541309,0.0716927973382896,0.0455239678491316,0.1420687775652751,0.0000000000000000,0.2417200033382098,0.0524420958026338,0.0019358291541309,0.0716927973382896,-0.0455239678491316,0.1420687775652751,0.0000000000000000,0.0524420958026338,0.2417200033382098,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0599992976854208,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.3858070096404518,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0599992976854208,0.3858070096404518,0.0000000000000000,-0.0000000000000000,0.0000000000000000,0.1065013559313907,0.1065013559313907,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.1065013559313907,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.1065013559313907,0.0000000000000000,0.0000000000000000,0.1359228235778018,-0.0011430477459907,-0.0192919647703918,-0.0150726194915292,0.0000000000000000,-0.0012527366048695,0.0002314796635360,-0.0011430477459906,-0.0952400040728509,-0.0621419541292237,-0.0485508884241462,0.0000000000000000,0.0135499675178389,-0.0141923898249493,-0.0192919647703918,-0.0621419541292237,-0.0782818474951750,0.0283405565751605,0.0000000000000000,0.0313151896932430,-0.0039145372579601,-0.0150726194915292,-0.0485508884241462,0.0283405565751605,-0.0924137079306580,0.0000000000000000,0.0244662451009505,-0.0162422252766831,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.1145559018709632,0.0000000000000000,0.0000000000000000,-0.0012527366048695,0.0135499675178390,0.0313151896932430,0.0244662451009505,0.0000000000000000,0.4208492761946447,0.0206525319754159,0.0002314796635360,-0.0141923898249491,-0.0039145372579601,-0.0162422252766831,0.0000000000000000,0.0206525319754159,0.0411253530731158,0.1016996853369585,-0.0008552470661841,0.0257839374356668,-0.0112775810503584,0.0000000000000000,0.0006097756403931,-0.0013738968295203,-0.0008552470661840,-0.0712601327043167,0.0830534513445227,-0.0363265708112602,0.0000000000000000,-0.0187793286172432,0.0182986577132198,0.0257839374356668,0.0830534513445227,-0.0585717618888200,-0.0378774866283914,0.0000000000000000,-0.0051309171701028,-0.0314903765392903,-0.0112775810503584,-0.0363265708112602,-0.0378774866283914,-0.0691454515877523,0.0000000000000000,-0.0241270010815564,0.0302803475639016,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0000000000000000,-0.0857126041609919,0.0000000000000000,0.0000000000000000,0.0006097756403931,-0.0187793286172431,-0.0051309171701028,-0.0241270010815564,0.0000000000000000,-0.0809242445997505,0.0154525630650178,-0.0013738968295203,0.0182986577132198,-0.0314903765392903,0.0302803475639016,0.0000000000000000,0.0154525630650178,0.4265812372482036,
]


def run(mm):
    tester = psr.PyTester("Testing Building of the DF coefficients")
    wf=make_wf()
    mm.load_supermodule("pulsar_libint")
    mm.load_supermodule("pulsar_scf")
    mm.change_option("PSR_3C2E","DF_INTS_KEY","LIBINT_3C2E")
    mm.change_option("PSR_Metric","METRIC_INTS_KEY","LIBINT_Metric")
    bs=wf.system.get_basis_set("PRIMARY")
    ints=mm.get_module("PSR_DFCoef",0).calculate("???",0,wf,bs,bs,bs)[0]
    dims=ints.sizes()
    flat_ints=[]
    for i in range(dims[0]):
        for j in range(dims[1]):
            for k in range(dims[2]):
                flat_ints.append(ints.get_value([i,j,k]))
    tester.test_double_vector("d^Q_{mn}",flat_ints,corr_ints)
    tester.print_results()
    return tester.nfailed()

def run_test():
    with psr.ModuleAdministrator() as mm:
        return run(mm)
