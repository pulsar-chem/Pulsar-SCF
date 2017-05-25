#include "pulsar_scf/builders/T.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <memory>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

const string T_opt="T_INTS_KEY";

MatrixBuilder::HashType TElectronic::my_hash_(const string &,
                                               unsigned int deriv,
                                               const Wavefunction &wfn,
                                               const BasisSet &bs1,
                                               const BasisSet &bs2)
{
    auto T_key=options().get<string>(T_opt);
    auto TInts=create_child_from_option<TwoCenterIntegral>(T_opt);
    return T_key+TInts->my_hash(deriv,wfn,bs1,bs2);
}

ReturnType TElectronic:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value was not in cache and FORCE_CACHE=true");
    auto TInts=create_child_from_option<TwoCenterIntegral>(T_opt);
    return matrix_builder_kernel(deriv,wfn,bs1,bs2,*TInts);
}

}//End namespace
