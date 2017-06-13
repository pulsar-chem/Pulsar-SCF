#include "pulsar_scf/builders/D.hpp"
#include "pulsar_scf/builders/Builders.hpp"
#include "pulsar_scf/SCF.hpp"
#include "pulsar_scf/helpers/SchwarzScreen.hpp"
#include "pulsar_scf/ga_builders/Builders.hpp"
using pulsar::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs insert_supermodule(void){
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<pulsar_scf::CoreDensity>("CoreDensity");
    cf.add_cpp_creator<pulsar_scf::G>("G");
    cf.add_cpp_creator<pulsar_scf::F>("F");
    cf.add_cpp_creator<pulsar_scf::HCore>("HCore");
    cf.add_cpp_creator<pulsar_scf::JK>("JK");
    cf.add_cpp_creator<pulsar_scf::DFJK>("DFJK");
    cf.add_cpp_creator<pulsar_scf::Overlap>("Overlap");
    cf.add_cpp_creator<pulsar_scf::Metric>("Metric");
    cf.add_cpp_creator<pulsar_scf::DFInts>("DFInts");
    cf.add_cpp_creator<pulsar_scf::DFCoef>("DFCoef");
    cf.add_cpp_creator<pulsar_scf::TElectronic>("TElectronic");
    cf.add_cpp_creator<pulsar_scf::SCF>("SCF");
    cf.add_cpp_creator<pulsar_scf::SchwarzMetric>("SchwarzScreen");
    cf.add_cpp_creator<pulsar_scf::NuclearElectronic>("NuclearElectronic");
    cf.add_cpp_creator<pulsar_scf::Orthogonalizer>("Orthogonalizer");
    cf.add_cpp_creator<pulsar_scf::GAJK>("GAJK");
    cf.add_cpp_creator<pulsar_scf::GAT>("GAT");
    cf.add_cpp_creator<pulsar_scf::GAV>("GAV");
    cf.add_cpp_creator<pulsar_scf::GAS>("GAS");
    cf.add_cpp_creator<pulsar_scf::GAX>("GAX");
    cf.add_cpp_creator<pulsar_scf::GAH>("GAH");
    cf.add_cpp_creator<pulsar_scf::GAG>("GAG");
    cf.add_cpp_creator<pulsar_scf::GAF>("GAF");
    return cf;
}

void finalize_supermodule(){

}

}
