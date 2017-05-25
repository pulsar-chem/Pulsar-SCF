#include "pulsar_scf/builders/S.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <memory>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

const string s_opt="S_INTS_KEY";

MatrixBuilder::HashType Overlap::my_hash_(const string &,
                                           unsigned int deriv,
                                           const Wavefunction &wfn,
                                           const BasisSet &bs1,
                                           const BasisSet &bs2)
{
    auto s_key=options().get<string>(s_opt);
    auto SInts=create_child_from_option<TwoCenterIntegral>(s_opt);
    return s_key+SInts->my_hash(deriv,wfn,bs1,bs2);
}

ReturnType Overlap:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value was not in cache and FORCE_CACHE=true");
    auto SInts=create_child_from_option<TwoCenterIntegral>(s_opt);
    return matrix_builder_kernel(deriv,wfn,bs1,bs2,*SInts);
}

}//End namespace
