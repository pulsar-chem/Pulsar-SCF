#include "pulsar_scf/builders/G.hpp"
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

const string JK_opt="JK_KEY";

MatrixBuilder::HashType G::my_hash_(const string &key,
                                    unsigned int deriv,
                                    const Wavefunction &wfn,
                                    const BasisSet &bs1,
                                    const BasisSet &bs2)
{
    const string JK_key=options().get<string>(JK_opt);
    auto JK=create_child_from_option<MatrixBuilder>(JK_opt);
    return JK_key+JK->my_hash(key,deriv,wfn,bs1,bs2);
}
ReturnType G::calculate_(const string &,
                         unsigned int deriv,
                         const Wavefunction & wfn,
                         const BasisSet &bs1,
                         const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");
    const auto JK=create_child_from_option<MatrixBuilder>(JK_opt)
                    ->calculate("",deriv,wfn,bs1,bs2);
    const auto J=*convert_to_eigen(*JK[0]);
    const auto K=*convert_to_eigen(*JK[1]);
    const auto G=J+0.25*K;
    return {make_shared<EigenMatrixImpl>(move(G))};
}

}//End namespace
