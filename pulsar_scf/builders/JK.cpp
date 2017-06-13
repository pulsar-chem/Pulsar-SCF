#include "pulsar_scf/builders/Builders.hpp"
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

//Implements the hash for a direct JK builder
MatrixBuilder::HashType JK::my_hash_(const string &,
                                     unsigned int deriv,
                                     const Wavefunction &wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    //A direct JK builder depends on three things

    //The instance that makes the ERIs
    auto eri_key=options().get<string>(ERI_opt);

    //The density
    auto D_hash=bphash::hash_to_string(wfn.opdm->my_hash());

    //And the hash of the module that makes the ERIs
    auto eri=create_child_from_option<FourCenterIntegral>(ERI_opt);
    return eri_key+D_hash+eri->my_hash(deriv,wfn,bs1,bs2,bs1,bs2);
}

//Implements a direct JK builder
ReturnType JK::calculate_(const string &,
                          unsigned int deriv,
                          const Wavefunction & wfn,
                          const BasisSet &bs1,
                          const BasisSet &bs2)
{
    //Hack to ensure we are using the cache
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache, but FORCE_CACHE=true");
    //TODO: code up non-symmetric basis sets
    if(bs1!=bs2)
        throw PulsarException("Direct JK build does not support differt bs's");

    //Get the densit matrix
    auto D=*convert_to_eigen(*wfn.opdm->get(Irrep::A,Spin::alpha));

    //Function that will make ERIs
    const auto Ints=create_child_from_option<FourCenterIntegral>(ERI_opt);
    Ints->initialize(deriv,wfn,bs1,bs2,bs1,bs2);

    //Allocate memory for J and K
    const size_t nbf=bs1.n_functions();
    matrix_type J=matrix_type::Zero(nbf,nbf),
                K=matrix_type::Zero(nbf,nbf);

    //Metric for assessing the importance of a shell pair's contribution
    //to a quartet
    auto schwarz_metric=create_child<MatrixBuilder>("PSR_Sieve");
    const auto temp_metric=schwarz_metric->calculate("",deriv,wfn,bs1,bs2)[0];
    const auto metric=*convert_to_eigen(*temp_metric);
    //The object that actually does the sieving
    SchwarzScreen sieve(metric,D,bs1);

    //An iterator that iterates over the shell quartets
    ShellQuartetItr quarts(bs1);
    while(quarts)//Loop over shell quartets
    {
        //Our current quartet
        const auto& shell=*quarts;
        if(!sieve.is_good(shell))//Is the quartet screened out
        {
            ++quarts;
            continue;
        }

        //Get the number of times this quartet should be counted
        const double total_deg=quarts.degeneracy();
        //The actual quartet
        const double* buffer=
                   Ints->calculate(shell[0],shell[1],shell[2],shell[3]);

        //Digest the quartet
        size_t counter=0;
        for(const auto& bf_quart:quarts)//Loop over basis funcitons in quartet
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

    //Symmetrize J and K
    matrix_type J_final=0.5*(J+J.transpose());
    matrix_type K_final=0.5*(K+K.transpose());

    return {make_shared<EigenMatrixImpl>(move(J_final)),
            make_shared<EigenMatrixImpl>(move(K_final))};
}



//Implements the hash for a core DF-JK build
MatrixBuilder::HashType DFJK::my_hash_(const string& key,
                                 unsigned int deriv,
                                 const Wavefunction& wfn,
                                 const BasisSet& bs1,
                                 const BasisSet& bs2)
{
    //The hash of the JK build depends on three things:

    //The instance that builds the DF coefs
    auto coef_key=options().get<string>(Coef_opt);
    //The MO coefs
    auto C_hash=bphash::hash_to_string(wfn.cmat->my_hash());
    //The df basis (factors into hash of DF coef module)
    auto df_bs_key=options().get<string>(DF_bs_opt);
    auto dfbs=wfn.system->get_basis_set(df_bs_key);
    //The hash of the module that builds the DF coef
    auto ds=create_child_from_option<Rank3Builder>(Coef_opt);
    return coef_key+C_hash+ds->my_hash(key,deriv,wfn,dfbs,bs1,bs2);
}


//Implementation of a JK build that uses a core DF build
ReturnType DFJK::calculate_(const std::string & key,
                 unsigned int deriv,
                 const Wavefunction & wfn,
                 const BasisSet & bs1,
                 const BasisSet & bs2)
{
    //TODO: This is a core algorithm need a direct algorithm

    //Hack to force cacheing
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache, but FORCE_CACHE=true");
    //TODO: Handle different bs1 and bs2
    if(bs1!=bs2)
        throw PulsarException("DFJK does not support different basis sets");

    //This is the key for the fitting basis set
    const auto& dfbs_name=options().get<string>(DF_bs_opt);

    //Get the fitting basis set from the system
    const auto& dfbs=wfn.system->get_basis_set(dfbs_name);

    //Function that will make our DF coefs
    auto coef_generator=create_child_from_option<Rank3Builder>(Coef_opt);

    //The DF coefs D_{Qls}=[L^-1]_{QP}(P|ls) Q,P run over density fitting basis
    //functions and l,s over atomic orbitals
    const auto temp_coefs=
            coef_generator->calculate(key,deriv,wfn,dfbs,bs1,bs2)[0];
    const auto d_Qls=*convert_to_eigen(*temp_coefs);

    //The MO coefficients C_li l runs over atomic orbitals i runs over occupied
    //molecular orbitals
    const auto C=*convert_to_eigen(*wfn.cmat->get(Irrep::A,Spin::alpha));

    const size_t nbf1=bs1.n_functions(),nbf2=bs2.n_functions(),nocc=C.cols();

    //MO coefs in tensor API
    tensor_matrix C_temp(C.data(),nbf1,nocc);

    //Final resting places for J and K
    auto J=make_shared<matrix_type>(matrix_type::Zero(nbf1,nbf2)),
         K=make_shared<matrix_type>(matrix_type::Zero(nbf1,nbf2));

    //Apply tensor API to J and K
    Eigen::TensorMap<Eigen::Tensor<double,2>> J_temp(J->data(),nbf1,nbf2),
                                              K_temp(K->data(),nbf1,nbf2);
    using idx_t=std::pair<int,int>;

    //Common intermediate: d_iQS=C_li(Q|ls)
    //Contract over index 0 of tensor 0 and index 1 of tensor 1
    const std::array<idx_t,1> idx{std::make_pair(0,1)};
    tensor_type d_iQs=C_temp.contract(d_Qls,idx);

    //J intermediate: d_Q=C_si(Q|is)
    //we actually do this as d_Q=(Q|is)C_si (Eigen's in column major)
    //Contract over index 0 and 2 of tensor 0 and indices 1 and 0 of tensor 1
    const std::array<idx_t,2> idx2{std::make_pair(0,1),std::make_pair(2,0)};
    Eigen::Tensor<double,1> d_Q=d_iQs.contract(C_temp,idx2);

    //J_mn=d_Q(mn|Q), actually do this as d_Q(Q|mn) to avoid transposing
    //contract over index 0 of tensor 0 and index 0 of tensor 1
    const std::array<idx_t,1> idx3{std::make_pair(0,0)};
    J_temp+=d_Q.contract(d_Qls,idx3);

    //K_mn=(mi|Q)(Q|in), actually do this as (Q|im)(Q|in)
    //contract over index 0 and 1 of both tensors
    const std::array<idx_t,2> idx4{std::make_pair(0,0),std::make_pair(1,1)};
    K_temp-=d_iQs.contract(d_iQs,idx4);

    //TODO: worry about coefficients so not hard-coded to RHF
    (*J)*=2.0;
    (*K)*=4.0;

    return {make_shared<EigenMatrixImpl>(J),
            make_shared<EigenMatrixImpl>(K)};

}



}//End namespace
