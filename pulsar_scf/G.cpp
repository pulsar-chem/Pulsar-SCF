#include "pulsar_scf/G.hpp"
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=OneElectronMatrix::ReturnType;
namespace pulsar_scf {

ReturnType G::calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    const auto JK=create_child_from_option<OneElectronMatrix>("JK_KEY")
                    ->calculate("",deriv,wfn,bs1,bs2);
    const auto J=*convert_to_eigen(*JK[0]);
    const auto K=*convert_to_eigen(*JK[1]);
    const auto G=J+0.25*K;
    return {std::make_shared<EigenMatrixImpl>(std::move(G))};
}

}//End namespace
