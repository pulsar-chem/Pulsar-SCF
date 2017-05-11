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

void DFJK::make_coefs(unsigned int deriv,
                      const Wavefunction& wfn,
                      const BasisSet& dfbs,
                      const BasisSet& bs1,
                      const BasisSet& bs2)
{
    const size_t nbf1=bs1.n_functions(),nbf2=bs2.n_functions(),ndf=dfbs.n_functions();
    const auto metric_ints=
            create_child_from_option<MatrixBuilder>("METRIC_KEY");
    const auto df_ints=
            create_child_from_option<Rank3Builder>("DF_INTS_KEY");
    auto Jmetric =
            *convert_to_eigen(*metric_ints->calculate("",deriv,wfn,bs1,bs2)[0]);
    matrix_type Linv_temp=
          Jmetric.llt().matrixL().solve(matrix_type::Identity(nbf1,nbf2));
    tensor_type Pls=
          *convert_to_eigen(*df_ints->calculate("",deriv,wfn,dfbs,bs1,bs2)[0]);
    tensor_matrix Linv(Linv_temp.data(), ndf,ndf);
    std::array<Eigen::IndexPair<int>,1> idx({Eigen::IndexPair<int>(1,0)});
    d_Qls_=std::make_unique<Eigen::Tensor<double,3>>(Linv.contract(Pls,idx));
}

ReturnType DFJK::calculate_(const std::string & key,
                 unsigned int deriv,
                 const Wavefunction & wfn,
                 const BasisSet & bs1,
                 const BasisSet & bs2)
{
    const size_t nbf1=bs1.n_functions(),nbf2=bs2.n_functions();    
    const BasisSet& dfbs=wfn.system->get_basis_set("PRIMARY");
    const size_t ndf=dfbs.n_functions();
    if(!d_Qls_)make_coefs(deriv,wfn,dfbs,bs1,bs2);

    const auto C=*convert_to_eigen(*wfn.cmat->get(Irrep::A,Spin::alpha));
    const size_t nocc=C.cols();
    tensor_matrix C_temp(C.data(),nbf1,nocc);
    std::array<Eigen::IndexPair<int>,1> idx({Eigen::IndexPair<int>(1,0)});
    tensor_type d_Qis=d_Qls_->contract(C_temp,idx);
    for(size_t Q=0;Q<ndf;++Q)
        for(size_t i=0;i<nocc;++i)
            for(size_t l=0;l<nbf1;++l)
                std::cout<<d_Qis(Q,i,l)<<std::endl;
    matrix_type J=matrix_type::Zero(nbf1,nbf2),
                K=matrix_type::Zero(nbf1,nbf2);
    //tensor_matrix J_temp(J.data(),nbf1,nbf2),K_temp(K.data(),nbf1,nbf2);
    //J_temp=C_temp.contract()

    for(size_t mu=0;mu<nbf1;++mu)
        for(size_t nu=0;nu<nbf2;++nu)
            for(size_t lambda=0;lambda<nbf1;++lambda)
                for(size_t sigma=0;sigma<nbf2;++sigma)
                    for(size_t Q=0;Q<ndf;++Q)
                        for(size_t i=0;i<C.cols();++i)
                    {
                        J(mu,nu)+=C(lambda,i)*C.transpose()(i,sigma)*(*d_Qls_)(Q,mu,nu)*(*d_Qls_)(Q,lambda,sigma);
                        K(mu,nu)-=C(lambda,i)*C.transpose()(i,sigma)*(*d_Qls_)(Q,mu,sigma)*(*d_Qls_)(Q,lambda,nu);
                    }
    //std::cout<<J<<std::endl;
    //std::cout<<K<<std::endl;
    exit(0);

}



}//End namespace
