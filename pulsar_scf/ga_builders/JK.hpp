#pragma once
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include "pulsar_scf/ga_builders/GlobalArrays.hpp"

namespace pulsar_scf {
#ifdef ENABLE_GA
class GAJK: public pulsar::MatrixBuilder {
public:
    GAJK(ID_t id):MatrixBuilder(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);

    HashType my_hash_(const std::string & key,
                       unsigned int deriv,
                      const pulsar::Wavefunction &wfn,
                      const pulsar::BasisSet &bs1,
                      const pulsar::BasisSet &bs2);
};
#else
NO_GA_MATRIX_BUILDER(GAJK)
#endif


}//End namespace pulsar_scf
