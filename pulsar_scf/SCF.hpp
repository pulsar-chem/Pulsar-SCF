#pragma once
#include <pulsar/modulebase/EnergyMethod.hpp>

namespace pulsar_scf {

class SCF: public pulsar::EnergyMethod {
public:
    SCF(ID_t id):EnergyMethod(id){}
    pulsar::DerivReturnType deriv_(size_t deriv,
                          const pulsar::Wavefunction & wfn);
};

}
