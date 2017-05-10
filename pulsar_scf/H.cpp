#include "pulsar_scf/H.hpp"
#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=OneElectronMatrix::ReturnType;
namespace pulsar_scf {

ReturnType HCore:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto terms=options().get<std::vector<std::string>>("H_KEYS");
    matrix_type H=Eigen::MatrixXd::Zero(bs1.n_functions(),bs2.n_functions());
    for(const auto& ti: terms)
    {
        auto termi=create_child<OneElectronMatrix>(ti);
        H+=*convert_to_eigen(*termi->calculate("",deriv,wfn,bs1,bs2)[0]);
    }

    return {std::make_shared<EigenMatrixImpl>(std::move(H))};
}

}//End namespace
