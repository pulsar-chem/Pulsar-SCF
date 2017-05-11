#pragma once
#include <pulsar/modulebase/MatrixBuilder.hpp>

namespace pulsar_scf {

class HCore: public pulsar::MatrixBuilder {
public:
    HCore(ID_t id):MatrixBuilder(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
};

}
