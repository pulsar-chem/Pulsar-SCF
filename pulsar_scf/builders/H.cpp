#include "pulsar_scf/builders/H.hpp"
#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

MatrixBuilder::HashType HCore::my_hash_(const string & key,
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


ReturnType HCore:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");
    auto terms=options().get<vector<string>>("H_KEYS");
    matrix_type H=Eigen::MatrixXd::Zero(bs1.n_functions(),bs2.n_functions());
    for(const auto& ti: terms)
    {
        auto termi=create_child<MatrixBuilder>(ti);
        H+=*convert_to_eigen(*termi->calculate("",deriv,wfn,bs1,bs2)[0]);
    }

    return {make_shared<EigenMatrixImpl>(move(H))};
}

}//End namespace
