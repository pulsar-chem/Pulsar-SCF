#include "pulsar_scf/ga_builders/Builders.hpp"
#include "pulsar_scf/ga_builders/GATensor.hpp"

using namespace pulsar;
using namespace std;

namespace pulsar_scf {

MatrixBuilder::HashType GAH::my_hash_(const string & key,
                                               unsigned int deriv,
                                               const Wavefunction &wfn,
                                               const BasisSet &bs1,
                                               const BasisSet &bs2)
{
    string final_hash="";
    for(auto key : options().get<vector<string>>("H_KEYS"))
    {
        auto term=create_child<MatrixBuilder>(key);
        final_hash+=key+term->my_hash(key,deriv,wfn,bs1,bs2);
    }
    return final_hash;
}


MatrixBuilder::ReturnType GAH:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");
    auto terms=options().get<vector<string>>("H_KEYS");
    const size_t nbf=bs1.n_functions();
    GATensor H(std::array<size_t,2>({nbf,nbf}),0.0);
    for(const auto& ti: terms)
    {
        auto termi=create_child<MatrixBuilder>(ti);
        GAaccumulate(H,1.0,*convert_to_GA(*termi->calculate("",deriv,wfn,bs1,bs2)[0]));
    }

    return {make_shared<GATensorImpl<2>>(H)};
}

}//End namespace
