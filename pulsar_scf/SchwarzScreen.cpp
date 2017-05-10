#include "pulsar_scf/SchwarzScreen.hpp"
#include <pulsar/modulebase/TwoElectronIntegral.hpp>

using namespace pulsar;

using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=OneElectronMatrix::ReturnType;

namespace pulsar_scf {

ReturnType SchwarzScreen::calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet & bs1,
                                     const BasisSet & bs2)
{
    auto ERIInts=create_child_from_option<TwoElectronIntegral>("ERI_INTS_KEY");
    ERIInts->options().change("THRESHOLD",0.0);
    ERIInts->initialize(deriv,wfn,bs1,bs2,bs1,bs2);
    const bool is_symmetric=bs1==bs2;
    const size_t n1=bs1.n_shell(),n2=bs2.n_shell();
    matrix_type S(bs1.n_shell(),bs2.n_shell());
    for(size_t i=0;i<n1;++i)
    {
        const size_t nbf_i=bs1.shell(i).n_functions();
        for(size_t j=(is_symmetric? i: 0);j<n2;++j)
        {
            const double* eri=ERIInts->calculate(i,j,i,j);
            const size_t nbf_j=bs2.shell(j).n_functions();
            Eigen::Map<const matrix_type> buffer(eri, nbf_i*nbf_j, nbf_i*nbf_j);
            const double value = std::sqrt(buffer.norm());
            S(i,j)=value;
            if(is_symmetric)S(j,i)=value;
        }
    }
    return {std::make_shared<EigenMatrixImpl>(std::move(S))};

}

}
