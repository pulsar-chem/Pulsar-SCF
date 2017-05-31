#include "pulsar_scf/ga_builders/Builders.hpp"
#include "pulsar_scf/ga_builders/GATensor.hpp"
#include <memory>

using namespace pulsar;
using namespace std;

namespace pulsar_scf {

const string H_opt="H_KEY";
const string G_opt="G_KEY";

MatrixBuilder::HashType GAF::my_hash_(const string &key,
                                    unsigned int deriv,
                                    const Wavefunction &wfn,
                                    const BasisSet &bs1,
                                    const BasisSet &bs2)
{
    const string H_key=options().get<string>(H_opt);
    const string G_key=options().get<string>(G_opt);
    auto H=create_child_from_option<MatrixBuilder>(H_opt);
    auto G=create_child_from_option<MatrixBuilder>(G_opt);
    return H_key+G_key+H->my_hash(key,deriv,wfn,bs1,bs2)+
            G->my_hash(key,deriv,wfn,bs1,bs2);
}

MatrixBuilder::ReturnType GAF::calculate_(const string &,
                         unsigned int deriv,
                         const Wavefunction & wfn,
                         const BasisSet &bs1,
                         const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");
    const auto H=*convert_to_GA(*create_child_from_option<MatrixBuilder>(H_opt)
                    ->calculate("",deriv,wfn,bs1,bs2)[0]);
    const auto G=*convert_to_GA(*create_child_from_option<MatrixBuilder>(G_opt)
                                   ->calculate("",deriv,wfn,bs1,bs2)[0]);
    GATensor F=GAAdd(1.0,H,1.0,G);
    return {make_shared<GATensorImpl<2>>(F)};
}

}//End namespace
