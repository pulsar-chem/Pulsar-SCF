#include "pulsar_scf/builders/Builders.hpp"
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;

namespace pulsar_scf {

//The key for the instance that will build JK for us
const string JK_opt="JK_KEY";

//Implements the hash for a G builder
MatrixBuilder::HashType G::my_hash_(const string &key,
                                    unsigned int deriv,
                                    const Wavefunction &wfn,
                                    const BasisSet &bs1,
                                    const BasisSet &bs2)
{
    //A G builder only depends on the JK builder under it
    //TODO: verify that ROHF/UHF don't change this
    const string JK_key=options().get<string>(JK_opt);
    auto JK=create_child_from_option<MatrixBuilder>(JK_opt);
    return JK_key+JK->my_hash(key,deriv,wfn,bs1,bs2);
}

//Implements a G build
ReturnType G::calculate_(const string &,
                         unsigned int deriv,
                         const Wavefunction & wfn,
                         const BasisSet &bs1,
                         const BasisSet &bs2)
{
    //Hack to verify cacheing
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");

    //The JK builder
    auto JK_builder=create_child_from_option<MatrixBuilder>(JK_opt);
    const auto JK=JK_builder->calculate("",deriv,wfn,bs1,bs2);

    //By convention J is matrix 0 and K is matrix 1
    //TODO: In Pulsar-Core make a JK module type?
    const auto J=*convert_to_eigen(*JK[0]);
    const auto K=*convert_to_eigen(*JK[1]);

    //Build G
    const auto G=J+0.25*K;

    return {make_shared<EigenMatrixImpl>(move(G))};
}

}//End namespace
