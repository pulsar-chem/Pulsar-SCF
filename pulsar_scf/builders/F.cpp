#include "pulsar_scf/builders/F.hpp"
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

ReturnType F::calculate_(const std::string &,
                         unsigned int deriv,
                         const Wavefunction & wfn,
                         const BasisSet &bs1,
                         const BasisSet &bs2)
{
    const auto H=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("H_KEY")
                    ->calculate("",deriv,wfn,bs1,bs2)[0]);
    const auto G=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("G_KEY")
                                   ->calculate("",deriv,wfn,bs1,bs2)[0]);
    const auto F=H+G;
    return {std::make_shared<EigenMatrixImpl>(std::move(F))};
}

}//End namespace
