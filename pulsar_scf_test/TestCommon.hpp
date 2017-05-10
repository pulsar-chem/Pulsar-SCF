#pragma once
#include<pulsar/datastore/Wavefunction.hpp>

struct ptr_wrapper{
    const double* ptr_;
    size_t n_;
    const double* begin()const{return ptr_;}
    const double* end()const{return ptr_+n_;}
};

inline pulsar::Wavefunction make_wf()
{
pulsar::AtomSetUniverse MyU;
const double angstrom_to_bohr = 1 / 0.52917721092;
const double Oy=-0.07579*angstrom_to_bohr;
const double Hx=0.86681*angstrom_to_bohr;
const double Hy=0.60144*angstrom_to_bohr;
auto cg = pulsar::ShellType::CartesianGaussian;
const pulsar::BasisShellInfo O1s(cg,0,3,1,
                         {130.709320000, 23.808861000, 6.443608300},
                         {0.15432897, 0.53532814, 0.44463454});
const pulsar::BasisShellInfo O2s(cg,0,3,1,
                         {5.033151300, 1.169596100, 0.380389000},
                         {-0.09996723, 0.39951283, 0.70011547});
const pulsar::BasisShellInfo O2p(cg,1,3,1,
                         {5.033151300, 1.169596100, 0.380389000},
                         {0.15591627, 0.60768372, 0.39195739});
const pulsar::BasisShellInfo H1s(cg,0,3,1,
                         {3.425250910, 0.623913730, 0.168855400},
                         {0.15432897, 0.53532814, 0.44463454});
pulsar::Atom O=pulsar::create_atom({0.00000,Oy,0.00000},8);
O.basis_sets["PRIMARY"].shells=std::vector<pulsar::BasisShellInfo>({O1s,O2s,O2p});
pulsar::Atom H1=pulsar::create_atom({Hx,Hy,0.00000},1);
H1.basis_sets["PRIMARY"].shells.push_back(H1s);
pulsar::Atom H2=pulsar::create_atom({-Hx,Hy,0.00000},1);
H2.basis_sets["PRIMARY"].shells.push_back(H1s);
MyU.insert(O);
MyU.insert(H1);
MyU.insert(H2);
pulsar::Wavefunction wf;
wf.system=std::make_shared<pulsar::System>(MyU,true);
return wf;
}
