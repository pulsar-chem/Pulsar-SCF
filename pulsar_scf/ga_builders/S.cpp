#include "pulsar_scf/ga_builders/Builders.hpp"
#include "pulsar_scf/ga_builders/GlobalArrays.hpp"
#include "pulsar_scf/ga_builders/GATensor.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>

using namespace pulsar;
using namespace std;

namespace pulsar_scf {

const string T_opt="S_INTS_KEY";

MatrixBuilder::HashType GAS::my_hash_(const string &,
                                               unsigned int deriv,
                                               const Wavefunction &wfn,
                                               const BasisSet &bs1,
                                               const BasisSet &bs2)
{
    auto T_key=options().get<string>(T_opt);
    auto TInts=create_child_from_option<TwoCenterIntegral>(T_opt);
    return T_key+TInts->my_hash(deriv,wfn,bs1,bs2);
}

MatrixBuilder::ReturnType GAS:: calculate_(const string &,
                                           unsigned int deriv,
                                           const Wavefunction & wfn,
                                           const BasisSet &bs1,
                                           const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value was not in cache and FORCE_CACHE=true");
    auto TInts=create_child_from_option<TwoCenterIntegral>(T_opt);
    TInts->initialize(deriv,wfn,bs1,bs2);
    auto rv = make_shared<GATensorImpl<2>>(fill_symmetric(TInts,bs1));
    return {rv};
}

}//End namespace
