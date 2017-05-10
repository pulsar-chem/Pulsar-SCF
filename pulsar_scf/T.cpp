#include "pulsar_scf/T.hpp"
#include "pulsar_scf/MatrixFillFxns.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <memory>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=OneElectronMatrix::ReturnType;
namespace pulsar_scf {

ReturnType TElectronic:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto TInts=create_child_from_option<TwoCenterIntegral>("T_INTS_KEY");
    TInts->initialize(deriv,wfn,bs1,bs2);
    const bool is_symmetric= bs1==bs2;
    auto T(std::move(is_symmetric ?
                       detail_::fill_symmetric<matrix_type>(TInts,bs1) :
                       detail_::fill_asymmetric<matrix_type>(TInts,bs1,bs2)
           ));
    return {std::make_shared<EigenMatrixImpl>(std::move(T))};
}

}//End namespace
