#pragma once
#include <pulsar/modulebase/OneElectronMatrix.hpp>

namespace pulsar_scf {

class Orthogonalizer: public pulsar::OneElectronMatrix {
public:
    Orthogonalizer(ID_t id):OneElectronMatrix(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
};

}
