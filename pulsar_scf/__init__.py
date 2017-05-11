#from .pulsar_scf import *
from .modinfo import *

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
    mm.load_module("pulsar_scf","SCF","PSR_SCF")
    mm.load_module("pulsar_scf","SchwarzScreen","PSR_Sieve")
