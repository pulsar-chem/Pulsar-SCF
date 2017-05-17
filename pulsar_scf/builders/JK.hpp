#pragma once
#include <pulsar/modulebase/MatrixBuilder.hpp>
namespace pulsar_scf {

class JK: public pulsar::MatrixBuilder {
public:
    JK(ID_t id):MatrixBuilder(id){}
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

class DFJK: public pulsar::MatrixBuilder {
public:
    DFJK(ID_t id):MatrixBuilder(id){}
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


}
