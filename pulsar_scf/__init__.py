from .modinfo import *
try: #See if we need to initialize stuff
    from .pulsar_scf import *
    GA_initialize()

except:
    pass

def initialize(mm):
    mm.load_module("pulsar_scf","DFInts","PSR_3C2E")
    mm.load_module("pulsar_scf","Metric","PSR_Metric")
    mm.load_module("pulsar_scf","TElectronic","PSR_T")
    mm.load_module("pulsar_scf","NuclearElectronic","PSR_V")
    mm.load_module("pulsar_scf","F","PSR_F")
    mm.load_module("pulsar_scf","G","PSR_G")
    mm.load_module("pulsar_scf","HCore","PSR_H")
    mm.load_module("pulsar_scf","Overlap","PSR_S")
    mm.load_module("pulsar_scf","Orthogonalizer","PSR_X")
    mm.load_module("pulsar_scf","CoreDensity","PSR_DCore")
    mm.load_module("pulsar_scf","JK","PSR_JK")
    mm.load_module("pulsar_scf","DFJK","PSR_DFJK")
    mm.load_module("pulsar_scf","DFCoef","PSR_DFCoef")
    mm.load_module("pulsar_scf","SCF","PSR_SCF")
    mm.load_module("pulsar_scf","SchwarzScreen","PSR_Sieve")
    mm.load_module("pulsar_scf","GAJK","PSR_GA_JK")
    mm.load_module("pulsar_scf","GAT","PSR_GAT")
    mm.load_module("pulsar_scf","GAV","PSR_GAV")
    mm.load_module("pulsar_scf","GAS","PSR_GAS")
    mm.load_module("pulsar_scf","GAX","PSR_GAX")
    mm.load_module("pulsar_scf","GAH","PSR_GAH")
    mm.load_module("pulsar_scf","GAG","PSR_GAG")
    mm.load_module("pulsar_scf","GAF","PSR_GAF")
