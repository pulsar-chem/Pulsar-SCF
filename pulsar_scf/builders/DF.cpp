#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <pulsar/modulebase/ThreeCenterIntegral.hpp>
#include "pulsar_scf/builders/Builders.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include "pulsar_scf/helpers/ShellTripleItr.hpp"

using namespace pulsar;
using namespace std;
//Type of an Eigen matrix
using matrix_type=EigenMatrixImpl::matrix_type;
//Type of a rank 3 Eigen tensor
using tensor_type=EigenTensorImpl<3>::tensor_type;
//Type of the raw data from an Eigen Matrix wrapped in the Eigen Tensor API
using tensor_matrix=Eigen::TensorMap<Eigen::Tensor<const double,2>>;


namespace pulsar_scf {

//Option keys we will nedd, in one common place

//Key for the module that makes (P|Q)
const string m_opt="METRIC_INTS_KEY";
//Key for the module that turns (P|Q) into a matrix
const string metric_opt="METRIC_KEY";
//Key for the module that makes (P|mn)
const string df_opt="DF_INTS_KEY";

//Implements the hash function for the metric matrix builder
MatrixBuilder::HashType Metric::my_hash_(const string & key,
                                         unsigned int deriv,
                                         const Wavefunction & wfn,
                                         const BasisSet &bs1,
                                         const BasisSet &bs2)
{
    //The Metric depends on three things:
    //Which module instance made the integrals
    auto m_key=options().get<string>(m_opt);
    //and the hash of that instance
    auto MInts=create_child_from_option<TwoCenterIntegral>(m_opt);
    return m_key+MInts->my_hash(deriv,wfn,bs1,bs2);
}

//Implements the metric matrix builder
MatrixBuilder::ReturnType Metric::calculate_(const string &,
                                             unsigned int deriv,
                                             const Wavefunction & wfn,
                                             const BasisSet &bs1,
                                             const BasisSet &bs2)
{
    //This is a function that can make shell pairs (P|Q)
    auto MInts=create_child_from_option<TwoCenterIntegral>(m_opt);

    //Matrix builder kernel is factored out code for copying integral batches
    //into a matrix
    return matrix_builder_kernel(deriv,wfn,bs1,bs2,*MInts);
}

//Implements the hash for the (P|mn) integrals
Rank3Builder::HashType DFInts::my_hash_(const string &,
                                        unsigned int deriv,
                                        const Wavefunction& wfn,
                                        const BasisSet &bs1,
                                        const BasisSet &bs2,
                                        const BasisSet &bs3)
{
    //The matrix of DF integrals, (P|mn) depends on only two things:
    //The module instance that made the integrals
    auto df_key=options().get<string>(df_opt);
    //and the hash of that instance
    auto df_ints=create_child_from_option<ThreeCenterIntegral>(df_opt);
    return df_key+df_ints->my_hash(deriv,wfn,bs1,bs2,bs3);
}

//Implements the building of the (P|mn) integral matrix
Rank3Builder::ReturnType DFInts::calculate_(const string &,
                                            unsigned int deriv,
                                            const pulsar::Wavefunction & wfn,
                                            const pulsar::BasisSet & bs1,
                                            const pulsar::BasisSet & bs2,
                                            const pulsar::BasisSet & bs3)
{
    //This is a hack to ensure that cacheing works
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Cache does not contain value and FORCE_CACHE=True");

    if(bs2!=bs3)
        throw PulsarException("Non-symmetric case not coded");

    //The function that actually makes our integrals
    auto DFInts=
            create_child_from_option<ThreeCenterIntegral>(df_opt);
    DFInts->initialize(deriv,wfn,bs1,bs2,bs3);

    //This iterator will loop over all possible shell triplets (P|mn) with
    //P taken from bs1 and m and n taken from bs2
    //TODO: worry about bs3
    ShellTripleItr shell_triple(bs1,bs2);
    const size_t nbf1=bs1.n_functions(),
                 nbf2=bs2.n_functions();

    //Allocates the tensor that we will return
    Eigen::Tensor<double,3> ints(nbf1,nbf2,nbf2);

    //Now fill said tensor, first loop over shell triplets
    while(shell_triple)
    {
        //Get our current index, (P|mn)
        const auto& idx=*shell_triple;
        //Compute the current shell triplet of integrals
        auto buffer=DFInts->calculate(idx[0],idx[1],idx[2]);
        //If it is screened out continue
        if(buffer==nullptr){
            ++shell_triple;
            continue;
        }

        //Now we loop over basis functions in the shell triplet in order to
        //copy elements over to the tensor, taking advantage of symmetry
        //Shell triplets come out in row-major form so we can just burn a
        //counter off rather than computing an actual index
        size_t counter=0;
        for(const auto& bf_triple:shell_triple)
        {
            const size_t P=bf_triple[0];
            const size_t m=bf_triple[1];
            const size_t n=bf_triple[2];
            ints(P,m,n)=ints(P,n,m)=buffer[counter++];
        }

        //Next shell triplet
        ++shell_triple;
    }

    return {make_shared<EigenTensorImpl<3>>(move(ints))};
}

//Implements the hash function for computing the DF coefficients
Rank3Builder::HashType DFCoef::my_hash_(const string & key,
                                        unsigned int deriv,
                                        const Wavefunction& wfn,
                                        const BasisSet &bs1,
                                        const BasisSet &bs2,
                                        const BasisSet &bs3)
{
    //The DF coefficients depend on four things:

    //The module instance used to form (P|mn)
    auto df_key=options().get<string>(df_opt);

    //The module instance used to form (P|Q)
    auto metric_key=options().get<string>(metric_opt);

    //The hash of the module instance that formed (P|mn)
    auto df_ints=create_child_from_option<Rank3Builder>(df_opt);
    auto df_part=df_key+df_ints->my_hash(key,deriv,wfn,bs1,bs2,bs3);

    //The hash of the module instance that formed (P|Q)
    auto metric=create_child_from_option<MatrixBuilder>(metric_opt);
    auto metric_part=metric_key+metric->my_hash(key,deriv,wfn,bs1,bs1);

    return df_part+metric_part;
}

Rank3Builder::ReturnType DFCoef::calculate_(const std::string&,
                                            unsigned int deriv,
                                            const Wavefunction& wfn,
                                            const pulsar::BasisSet& dfbs,
                                            const pulsar::BasisSet& bs1,
                                            const pulsar::BasisSet& bs2)
{
    //Hack to ensure caching is working
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Cache value not found and FORCE_CACHE=True");

    //Number of df basis functions
    const size_t ndf=dfbs.n_functions();
    //The function that will make (P|Q)
    const auto metric_ints=
            create_child_from_option<MatrixBuilder>(metric_opt);
    //The function that will make (P|mn)
    const auto df_ints=
            create_child_from_option<Rank3Builder>(df_opt);

    //Now we invert (P|Q), we do this via LL^T (Cholesky) decomposition
    //This is (P|Q)
    const auto Jmetric_temp=metric_ints->calculate("",deriv,wfn,dfbs,dfbs)[0];
    const auto Jmetric=*convert_to_eigen(*Jmetric_temp);

    //The part up to .matrixL() is L of the LL^T decomposition of  (P|Q)
    //L^-1 is X such that LX=1
    matrix_type Linv_temp=
            Jmetric.llt().matrixL().solve(matrix_type::Identity(ndf,ndf));
    //L^-1 in the tensor API
    tensor_matrix Linv(Linv_temp.data(), ndf,ndf);
    //These are our (P|mn) integrals
    const auto Pmn_temp=df_ints->calculate("",deriv,wfn,dfbs,bs1,bs2)[0];
    tensor_type Pmn=*convert_to_eigen(*Pmn_temp);

    //Eigen tensor needs an array of indices to contract. This says we are
    //contracting index 1 of tensor 0 with index 0 of tensor 1
    std::array<std::pair<int,int>,1> idx{std::make_pair(1,0)};

    //This is the contraction [L^-1]_{PQ}(Q|mn)
    auto coefs=Linv.contract(Pmn,idx);

    return {std::make_shared<EigenTensorImpl<3>>(move(coefs))};
}

}//End namespace
