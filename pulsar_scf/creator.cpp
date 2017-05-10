#include "pulsar_scf/D.hpp"
#include "pulsar_scf/DF.hpp"
#include "pulsar_scf/G.hpp"
#include "pulsar_scf/F.hpp"
#include "pulsar_scf/H.hpp"
#include "pulsar_scf/JK.hpp"
#include "pulsar_scf/S.hpp"
#include "pulsar_scf/SCF.hpp"
#include "pulsar_scf/SchwarzScreen.hpp"
#include "pulsar_scf/T.hpp"
#include "pulsar_scf/V.hpp"
#include "pulsar_scf/X.hpp"

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
    cf.add_cpp_creator<pulsar_scf::TElectronic>("TElectronic");
    cf.add_cpp_creator<pulsar_scf::SCF>("SCF");
    cf.add_cpp_creator<pulsar_scf::SchwarzScreen>("SchwarzScreen");
    cf.add_cpp_creator<pulsar_scf::NuclearElectronic>("NuclearElectronic");
    cf.add_cpp_creator<pulsar_scf::Orthogonalizer>("Orthogonalizer");
    return cf;
}

}
