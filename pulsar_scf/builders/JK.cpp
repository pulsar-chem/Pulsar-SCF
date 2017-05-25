#include "pulsar_scf/builders/JK.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/FourCenterIntegral.hpp>
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include <pulsar/modulebase/Rank3Builder.hpp>
#include <memory>
#include <pulsar/math/EigenImpl.hpp>
#include "pulsar_scf/helpers/ShellQuartetItr.hpp"
#include "pulsar_scf/helpers/SchwarzScreen.hpp"

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using tensor_type=EigenTensorImpl<3>::tensor_type;
using tensor_matrix=Eigen::TensorMap<Eigen::Tensor<const double,2>>;
using ReturnType=MatrixBuilder::ReturnType;

namespace pulsar_scf {

const string ERI_opt="ERI_KEY";
const string DF_bs_opt="FITTING_BASIS_KEY";
const string Coef_opt="FITTING_COEF_KEY";


MatrixBuilder::HashType JK::my_hash_(const string &,
                                     unsigned int deriv,
                                     const Wavefunction &wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto eri_key=options().get<string>(ERI_opt);
    auto eri=create_child_from_option<FourCenterIntegral>(ERI_opt);
    auto D_hash=bphash::hash_to_string(wfn.opdm->my_hash());
    return eri_key+D_hash+eri->my_hash(deriv,wfn,bs1,bs2,bs1,bs2);
}

ReturnType JK::calculate_(const string &,
                          unsigned int deriv,
                          const Wavefunction & wfn,
                          const BasisSet &bs1,
                          const BasisSet &bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache, but FORCE_CACHE=true");
    auto sharedD=convert_to_eigen(*wfn.opdm->get(Irrep::A,Spin::alpha));
    const auto& D=*sharedD;
    const auto Ints=create_child_from_option<FourCenterIntegral>(ERI_opt);
    Ints->initialize(deriv,wfn,bs1,bs2,bs1,bs2);
    const size_t nbf=bs1.n_functions();
    matrix_type J=matrix_type::Zero(nbf,nbf),K=matrix_type::Zero(nbf,nbf);
    auto schwarz_metric=create_child<MatrixBuilder>("PSR_Sieve");
    const auto& metric=
         *convert_to_eigen(*schwarz_metric->calculate("",deriv,wfn,bs1,bs2)[0]);
    SchwarzScreen sieve(metric,D,bs1);
    ShellQuartetItr quarts(bs1);
    while(quarts)
    {
        const auto& shell=*quarts;
        if(!sieve.is_good(shell))
        {
            ++quarts;
            continue;
        }
        const double total_deg=quarts.degeneracy();
        const double* buffer=
                   Ints->calculate(shell[0],shell[1],shell[2],shell[3]);

        size_t counter=0;
        for(const auto& bf_quart:quarts)
        {
            const double value=buffer[counter++]*total_deg;
            const size_t mu=bf_quart[0],
                         nu=bf_quart[1],
                         lambda=bf_quart[2],
                         sigma=bf_quart[3];
            J(mu,nu)+=D(lambda,sigma)*value;
            J(lambda,sigma)+=D(mu,nu)*value;
            K(mu,lambda)-=D(nu,sigma)*value;
            K(mu,sigma) -=D(nu,lambda)*value;
            K(nu,sigma) -=D(mu,lambda)*value;
            K(nu,lambda)-=D(mu,sigma)*value;
        }
        ++quarts;
    }

    matrix_type J_final=0.5*(J+J.transpose());
    matrix_type K_final=0.5*(K+K.transpose());

    return {make_shared<EigenMatrixImpl>(move(J_final)),
            make_shared<EigenMatrixImpl>(move(K_final))};
}

using idx_t=Eigen::IndexPair<int>;
template<size_t n> using idx_array=array<idx_t,n>;


MatrixBuilder::HashType DFJK::my_hash_(const string& key,
                                 unsigned int deriv,
                                 const Wavefunction& wfn,
                                 const BasisSet& bs1,
                                 const BasisSet& bs2)
{
    auto df_bs_key=options().get<string>(DF_bs_opt);
    auto dfbs=wfn.system->get_basis_set(df_bs_key);
    auto coef_key=options().get<string>(Coef_opt);
    auto C_hash=bphash::hash_to_string(wfn.cmat->my_hash());
    auto ds=create_child_from_option<Rank3Builder>(Coef_opt);
    return coef_key+C_hash+ds->my_hash(key,deriv,wfn,dfbs,bs1,bs2);
}

ReturnType DFJK::calculate_(const std::string & key,
                 unsigned int deriv,
                 const Wavefunction & wfn,
                 const BasisSet & bs1,
                 const BasisSet & bs2)
{
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache, but FORCE_CACHE=true");
    const auto& dfbs_name=options().get<string>(DF_bs_opt);
    const BasisSet& dfbs=wfn.system->get_basis_set(dfbs_name);
    auto coef_generator=create_child_from_option<Rank3Builder>(Coef_opt);
    const auto d_Qls=*convert_to_eigen(
                *coef_generator->calculate(key,deriv,wfn,dfbs,bs1,bs2)[0]
    );
    //TODO: Handle different C left and C right
    const auto C=*convert_to_eigen(*wfn.cmat->get(Irrep::A,Spin::alpha));

    const size_t nbf1=bs1.n_functions(),nbf2=bs2.n_functions(),nocc=C.cols();

    tensor_matrix C_temp(C.data(),nbf1,nocc);

    auto J=make_shared<matrix_type>(matrix_type::Zero(nbf1,nbf2)),
         K=make_shared<matrix_type>(matrix_type::Zero(nbf1,nbf2));

    Eigen::TensorMap<Eigen::Tensor<double,2>> J_temp(J->data(),nbf1,nbf2),
                                              K_temp(K->data(),nbf1,nbf2);

    //Common intermediate: d_iQS=C_li(Q|ls)
    tensor_type d_iQs=
            C_temp.contract(d_Qls,idx_array<1>({idx_t({0,1})}));

    //J intermediate: d_Q=C_si(Q|is)
    Eigen::Tensor<double,1> d_Q=
            d_iQs.contract(C_temp,idx_array<2>({idx_t({0,1}),idx_t({2,0})}));

    //J_mn=d_Q(mn|Q)
    J_temp+=d_Q.contract(d_Qls,idx_array<1>({idx_t({0,0})}));

    //K_mn=(mi|Q)(Q|in)
    K_temp-=d_iQs.contract(d_iQs,idx_array<2>({idx_t({0,0}),idx_t({1,1})}));

    //TODO: worry about coefficients so not hard-coded to RHF
    (*J)*=2.0;
    (*K)*=4.0;

    return {make_shared<EigenMatrixImpl>(J),
            make_shared<EigenMatrixImpl>(K)};

}



}//End namespace
