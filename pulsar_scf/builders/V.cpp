#include "pulsar_scf/builders/Builders.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <bphash/Hash.hpp>
#include <memory>

using namespace pulsar;
using namespace std;

using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

//The key for the electron-nuclear attraction integrals
const string V_opt="V_INTS_KEY";

//Implements hash for the electron-nulecus attraction integral matrix builder
MatrixBuilder::HashType NuclearElectronic::my_hash_(const string &,
                                               unsigned int deriv,
                                               const Wavefunction &wfn,
                                               const BasisSet &bs1,
                                               const BasisSet &bs2)
{
    //The hash for V only depends on the integrals
    auto V_key=options().get<string>(V_opt);
    auto VInts=create_child_from_option<TwoCenterIntegral>(V_opt);
    return V_key+VInts->my_hash(deriv,wfn,bs1,bs2);
}

//Implements the building of V
ReturnType NuclearElectronic:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    //Hack to ensure cacheing occurred
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value was not in cache and FORCE_CACHE=true");

    //Make a integral builder and defer to factored code.
    auto VInts=create_child_from_option<TwoCenterIntegral>(V_opt);
    return matrix_builder_kernel(deriv,wfn,bs1,bs2,*VInts);
}

}//End namespace
