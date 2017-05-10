#pragma once
#include <pulsar/modulebase/OneElectronMatrix.hpp>
namespace pulsar_scf {

class JK: public pulsar::OneElectronMatrix {
public:
    JK(ID_t id):OneElectronMatrix(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
};

class DFJK: public pulsar::OneElectronMatrix {
private:
    //std::unique_ptr<Eigen::MatrixXd> metric_;
    std::unique_ptr<Eigen::Tensor<double,3>> d_Qls_;
    void make_coefs(unsigned int deriv,
                          const pulsar::Wavefunction& wfn,
                          const pulsar::BasisSet& dfbs,
                          const pulsar::BasisSet& bs1,
                          const pulsar::BasisSet& bs2);
public:
    DFJK(ID_t id):OneElectronMatrix(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
};


}
