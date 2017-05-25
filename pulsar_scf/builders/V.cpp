#include "pulsar_scf/builders/V.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <bphash/Hash.hpp>
#include <memory>

using namespace pulsar;
using namespace std;

using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

const string V_opt="V_INTS_KEY";

MatrixBuilder::HashType NuclearElectronic::my_hash_(const string &,
                                               unsigned int deriv,
                                               const Wavefunction &wfn,
                                               const BasisSet &bs1,
                                               const BasisSet &bs2)
{
    auto V_key=options().get<string>(V_opt);
    auto VInts=create_child_from_option<TwoCenterIntegral>(V_opt);
    return V_key+VInts->my_hash(deriv,wfn,bs1,bs2);
}

ReturnType NuclearElectronic:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value was not in cache and FORCE_CACHE=true");
    auto VInts=create_child_from_option<TwoCenterIntegral>(V_opt);
    return matrix_builder_kernel(deriv,wfn,bs1,bs2,*VInts);
}

}//End namespace
