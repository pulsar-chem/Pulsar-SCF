#include "pulsar_scf/builders/D.hpp"
#include "pulsar_scf/builders/DF.hpp"
#include "pulsar_scf/builders/G.hpp"
#include "pulsar_scf/builders/F.hpp"
#include "pulsar_scf/builders/H.hpp"
#include "pulsar_scf/builders/JK.hpp"
#include "pulsar_scf/builders/S.hpp"
#include "pulsar_scf/SCF.hpp"
#include "pulsar_scf/helpers/SchwarzScreen.hpp"
#include "pulsar_scf/builders/T.hpp"
#include "pulsar_scf/builders/V.hpp"
#include "pulsar_scf/builders/X.hpp"
#include "pulsar_scf/ga_builders/JK.hpp"
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
    return cf;
}

void finalize_supermodule(){
    pulsar_scf::GA_finalize();
}

}
