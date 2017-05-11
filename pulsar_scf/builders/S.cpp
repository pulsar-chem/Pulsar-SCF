#include "pulsar_scf/builders/S.hpp"
#include "pulsar_scf/MatrixFillFxns.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <memory>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

ReturnType Overlap:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto SInts=create_child_from_option<TwoCenterIntegral>("S_INTS_KEY");
    SInts->initialize(deriv,wfn,bs1,bs2);
    const bool is_symmetric= bs1==bs2;
    auto S(std::move(is_symmetric ?
                       detail_::fill_symmetric<matrix_type>(SInts,bs1) :
                       detail_::fill_asymmetric<matrix_type>(SInts,bs1,bs2)
           ));
    return {std::make_shared<EigenMatrixImpl>(std::move(S))};
}

}//End namespace
