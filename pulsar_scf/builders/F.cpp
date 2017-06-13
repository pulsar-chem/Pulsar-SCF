#include "pulsar_scf/builders/Builders.hpp"
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;

namespace pulsar_scf {

//The key for the module that will make H
const string H_opt="H_KEY";
//The key for the module that will make G
const string G_opt="G_KEY";

//Implements the hash of the Fock builder instance
MatrixBuilder::HashType F::my_hash_(const string &key,
                                    unsigned int deriv,
                                    const Wavefunction &wfn,
                                    const BasisSet &bs1,
                                    const BasisSet &bs2)
{
    //A Fock builder only depends on two things

    //The instance that builds H
    const string H_key=options().get<string>(H_opt);
    auto H=create_child_from_option<MatrixBuilder>(H_opt);

    //The instance that builds G
    const string G_key=options().get<string>(G_opt);
    auto G=create_child_from_option<MatrixBuilder>(G_opt);

    return H_key+G_key+H->my_hash(key,deriv,wfn,bs1,bs2)+
            G->my_hash(key,deriv,wfn,bs1,bs2);
}

//Implements a Fock build
ReturnType F::calculate_(const string &,
                         unsigned int deriv,
                         const Wavefunction & wfn,
                         const BasisSet &bs1,
                         const BasisSet &bs2)
{

    //Hack to ensure cacheing works
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");

    //The H builder and the retrieval of H
    auto H_builder=create_child_from_option<MatrixBuilder>(H_opt);
    const auto H_temp=H_builder->calculate("",deriv,wfn,bs1,bs2)[0];
    const auto H=*convert_to_eigen(*H_temp);

    //The G builder and the retrieval of G
    auto G_builder=create_child_from_option<MatrixBuilder>(G_opt);
    const auto G_temp=G_builder->calculate("",deriv,wfn,bs1,bs2)[0];
    const auto G=*convert_to_eigen(*G_temp);

    //Build F
    const auto F=H+G;

    return {make_shared<EigenMatrixImpl>(move(F))};
}

}//End namespace
