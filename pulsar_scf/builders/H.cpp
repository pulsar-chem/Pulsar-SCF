#include "pulsar_scf/builders/Builders.hpp"
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

//Implements the has for the core Hamiltonian
MatrixBuilder::HashType HCore::my_hash_(const string & key,
                                               unsigned int deriv,
                                               const Wavefunction &wfn,
                                               const BasisSet &bs1,
                                               const BasisSet &bs2)
{
    //We loop over terms in the core Hamiltonian and add their hashes to ours
    //A core builder thus depends on whatever terms are in it
    string final_hash="";
    for(auto key : options().get<vector<string>>("H_KEYS"))
    {
        auto term=create_child<MatrixBuilder>(key);
        final_hash+=key+term->my_hash(key,deriv,wfn,bs1,bs2);
    }
    return final_hash;
}

//Implements the core builder
ReturnType HCore:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    //Hack to ensure that results are being cached
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");

    //This is a list of terms to include in the core Hamiltonian as set by the
    //user
    auto terms=options().get<vector<string>>("H_KEYS");

    //Allocate memory for H
    matrix_type H=Eigen::MatrixXd::Zero(bs1.n_functions(),bs2.n_functions());

    //Loop over terms in the list
    for(const auto& ti: terms)
    {
        auto termi=create_child<MatrixBuilder>(ti);
        const auto termi_temp=termi->calculate("",deriv,wfn,bs1,bs2)[0];
        H+=*convert_to_eigen(*termi_temp);
    }

    return {make_shared<EigenMatrixImpl>(move(H))};
}

}//End namespace
