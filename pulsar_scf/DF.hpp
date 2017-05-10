#pragma once
#include <pulsar/modulebase/OneElectronMatrix.hpp>
#include <pulsar/modulebase/Rank3Builder.hpp>

namespace pulsar_scf {

class Metric: public pulsar::OneElectronMatrix {
public:
    Metric(ID_t id):OneElectronMatrix(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
};

class DFInts: public pulsar::Rank3Builder {
public:
    DFInts(ID_t id):Rank3Builder(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2,
                          const pulsar::BasisSet & bs3);
};


}
