#include "pulsar_scf/helpers/SchwarzScreen.hpp"
#include <pulsar/modulebase/FourCenterIntegral.hpp>
#include "pulsar_scf/helpers/ShellPairItr.hpp"

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;

namespace pulsar_scf {

const string ERI_opt="ERI_INTS_KEY";

MatrixBuilder::HashType SchwarzMetric::my_hash_(const string & key,
                                                unsigned int deriv,
                                                const Wavefunction & wfn,
                                                const BasisSet & bs1,
                                                const BasisSet & bs2)
{
    auto eri_key=options().get<string>(ERI_opt);
    auto eris=create_child_from_option<FourCenterIntegral>(ERI_opt);
    return eri_key+eris->my_hash(deriv,wfn,bs1,bs2,bs1,bs2);
}


ReturnType SchwarzMetric::calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet & bs1,
                                     const BasisSet & bs2)
{
    auto ERIInts=create_child_from_option<FourCenterIntegral>(ERI_opt);
    ERIInts->options().change("THRESHOLD",0.0);
    ERIInts->initialize(deriv,wfn,bs1,bs2,bs1,bs2);
    const bool is_symmetric=bs1==bs2;
    if(!is_symmetric)
        throw PulsarException("Non-symmetric case is not coded yet");
    ShellPairItr shell_pair(bs1);
    matrix_type S(bs1.n_shell(false),bs1.n_shell(false));
    while(shell_pair)
    {
       const auto& idx=*shell_pair;
       const size_t i=idx[0],j=idx[1];
       const double* eri=ERIInts->calculate(i,j,i,j);
       const size_t nbf2=shell_pair.begin().size();
       Eigen::Map<const matrix_type> buffer(eri, nbf2, nbf2);
       const double value = std::sqrt(buffer.norm());
       S(i,j)=S(j,i)=value;
       ++shell_pair;
    }
    return {make_shared<EigenMatrixImpl>(move(S))};

}

SchwarzScreen::SchwarzScreen(const Eigen::MatrixXd& metric,
                             const Eigen::MatrixXd& density,
                             const BasisSet& bs,
                             double threshold):
    metric_(metric),
    density_(Eigen::MatrixXd::Zero(bs.n_shell(false),bs.n_shell(false))),
    threshold_(threshold)
{
    ShellPairItr shell_pairs(bs);
    while(shell_pairs)
    {
        const auto& pair=*shell_pairs;
        for(const auto& bf_pair : shell_pairs)
        {
            const double Dmn=density(bf_pair[0],bf_pair[1]);
            density_(pair[0],pair[1])+=Dmn*Dmn;
        }
        density_(pair[0],pair[1])=density_(pair[1],pair[0])=
                std::sqrt(density_(pair[0],pair[1]));
        ++shell_pairs;
    }
}

}
