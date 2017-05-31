#include "pulsar_scf/ga_builders/Builders.hpp"
#include "pulsar_scf/ga_builders/GlobalArrays.hpp"
#include "pulsar_scf/ga_builders/GATensor.hpp"

using namespace pulsar;
using namespace std;

namespace pulsar_scf {

const string S_opt="S_KEY";

MatrixBuilder::HashType GAX::my_hash_(const string &,
                                                  unsigned int deriv,
                                                  const Wavefunction &wfn,
                                                  const BasisSet &bs1,
                                                  const BasisSet &bs2)
{
    auto S_key=options().get<string>(S_opt);
    auto Sbuilder=create_child_from_option<MatrixBuilder>(S_opt);
    return S_key+Sbuilder->my_hash("",deriv,wfn,bs1,bs2);
}


///This is flat out stolen from libint2's hartree-fock++.cc routine
MatrixBuilder::ReturnType GAX:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");

    auto S=*convert_to_GA(*create_child_from_option<MatrixBuilder>(S_opt)
                             ->calculate("",deriv,wfn,bs1,bs2)[0]);
    auto eig_sys=EigenSolver(S);
    auto& values=eig_sys.first;
    auto& vectors=eig_sys.second;
    const bool symmetric=false;
    //TODO: worry about small eigenvalues
    GATensor sigma=vec2matrix(GApow(values,0.5));
    GATensor sigma_invsqrt = vec2matrix(GApow(values,-0.5));
    GATensor X=gemm(1.0,vectors,sigma_invsqrt);
    GATensor Xinv=gemm(1.0,vectors,sigma);
    // convert to symmetric, if needed
    if (symmetric) {
        //TODO: verify that this works (should be extracted as separate module)
        vectors.transpose();
        X = gemm(1.0,X,vectors);
        Xinv = gemm(1.0,Xinv,vectors);
    }
    return {make_shared<GATensorImpl<2>>(X),
            make_shared<GATensorImpl<2>>(Xinv)};
}

}//End namespace
