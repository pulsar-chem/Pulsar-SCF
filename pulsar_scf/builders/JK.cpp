#include "pulsar_scf/builders/JK.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/FourCenterIntegral.hpp>
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include <pulsar/modulebase/Rank3Builder.hpp>
#include <memory>
#include <pulsar/math/EigenImpl.hpp>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using tensor_type=EigenTensorImpl<3>::tensor_type;
using tensor_matrix=Eigen::TensorMap<Eigen::Tensor<const double,2>>;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {


ReturnType JK::calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{

    using ci=const size_t;
    const auto& D=*convert_to_eigen(*wfn.opdm->get(Irrep::A,Spin::alpha));
    const auto Ints=create_child_from_option<FourCenterIntegral>("ERI_KEY");
    Ints->initialize(deriv,wfn,bs1,bs1,bs1,bs1);
    ci nshells1=bs1.n_shell(), nbf=bs1.n_functions();

    matrix_type J=matrix_type::Zero(nbf,nbf),
                K=matrix_type::Zero(nbf,nbf);

    for(size_t shell_i=0; shell_i<nshells1;++shell_i)
    {
        ci nbf_i=bs1.shell(shell_i).n_functions();
        ci off_i=bs1.shell_start(shell_i);

        for(size_t shell_j=0; shell_j<=shell_i; ++shell_j)
        {
            ci nbf_j=bs1.shell(shell_j).n_functions();
            ci off_j=bs1.shell_start(shell_j);

            const bool ieqj(shell_i==shell_j);
            const double ij_deg(ieqj ? 1.0 : 2.0);

            for(size_t shell_k=0; shell_k<=shell_i;++shell_k)
            {
                ci nbf_k=bs1.shell(shell_k).n_functions();
                ci off_k=bs1.shell_start(shell_k);

                const bool ieqk(shell_i==shell_k);
                ci max_l( ieqk ? shell_j : shell_k);


                for(size_t shell_l=0; shell_l<=max_l;++shell_l)
                {
                    ci nbf_l=bs1.shell(shell_l).n_functions();
                    ci off_l=bs1.shell_start(shell_l);

                    const bool keql(shell_k==shell_l);
                    const bool jeql(shell_j==shell_l);
                    const double kl_deg(keql ? 1.0 : 2.0);
                    const double ij_kl_deg( ieqk && jeql? 1.0 : 2.0);
                    const double total_deg(ij_deg*kl_deg*ij_kl_deg);
                    //Buffer is nbfi by nbfj by nbfk by nbfl
                    const double* buffer=
                               Ints->calculate(shell_i,shell_j,shell_k,shell_l);
                    if(buffer == nullptr)continue;

                    for(size_t mu=0, counter=0;mu<nbf_i;++mu)
                    {
                        ci mu_large=mu+off_i;

                        for(size_t nu=0;nu<nbf_j;++nu)
                        {
                            ci nu_large=nu+off_j;

                            for(size_t lambda=0;lambda<nbf_k;++lambda)
                            {
                                ci lambda_large=lambda+off_k;

                                for(size_t sigma=0;sigma<nbf_l;++sigma,++counter)
                                {
                                    ci sigma_large=sigma+off_l;
                                    const double value=buffer[counter]*total_deg;

                                    J(mu_large,nu_large)+=D(lambda_large,sigma_large)*value;
                                    J(lambda_large,sigma_large)+=D(mu_large,nu_large)*value;
                                    K(mu_large,lambda_large)-=D(nu_large,sigma_large)*value;
                                    K(mu_large,sigma_large)-=D(nu_large,lambda_large)*value;
                                    K(nu_large,sigma_large)-=D(mu_large,lambda_large)*value;
                                    K(nu_large,lambda_large)-=D(mu_large,sigma_large)*value;
                                }
                            }
                        }
                    }

                }
            }
        }
    }

    matrix_type J_final=0.5*(J+J.transpose());
    matrix_type K_final=0.5*(K+K.transpose());

    return {std::make_shared<EigenMatrixImpl>(std::move(J_final)),
            std::make_shared<EigenMatrixImpl>(std::move(K_final))};
}

using idx_t=Eigen::IndexPair<int>;
template<size_t n> using idx_array=std::array<idx_t,n>;

ReturnType DFJK::calculate_(const std::string & key,
                 unsigned int deriv,
                 const Wavefunction & wfn,
                 const BasisSet & bs1,
                 const BasisSet & bs2)
{
    const auto& dfbs_name=options().get<std::string>("FITTING_BASIS_KEY");
    const BasisSet& dfbs=wfn.system->get_basis_set(dfbs_name);
    auto coef_generator=create_child_from_option<Rank3Builder>("FITTING_COEF_KEY");
    const auto d_Qls=*convert_to_eigen(
                *coef_generator->calculate("",deriv,wfn,dfbs,bs1,bs2)[0]
    );
    //TODO: Handle differently sized C left and C right
    const auto C=*convert_to_eigen(*wfn.cmat->get(Irrep::A,Spin::alpha));

    const size_t nbf1=bs1.n_functions(),nbf2=bs2.n_functions(),nocc=C.cols();

    tensor_matrix C_temp(C.data(),nbf1,nocc);

    auto J=std::make_shared<matrix_type>(matrix_type::Zero(nbf1,nbf2)),
         K=std::make_shared<matrix_type>(matrix_type::Zero(nbf1,nbf2));

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

    return {std::make_shared<EigenMatrixImpl>(J),
            std::make_shared<EigenMatrixImpl>(K)};

}



}//End namespace
