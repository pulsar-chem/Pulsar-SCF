#include "pulsar_scf/ga_builders/JK.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/FourCenterIntegral.hpp>
#include <pulsar/modulebase/Rank3Builder.hpp>
#include <memory>
#include "pulsar_scf/ga_builders/GATensor.hpp"
#include "pulsar_scf/helpers/ShellQuartetItr.hpp"
#include "pulsar_scf/helpers/SchwarzScreen.hpp"
#include <pulsar/math/EigenImpl.hpp>
#include <ga.h>

using namespace pulsar;
using namespace std;



using Idx_t = std::array<size_t,2>;

namespace pulsar_scf {

const string ERI_opt="ERI_KEY";
const string DF_bs_opt="FITTING_BASIS_KEY";
const string Coef_opt="FITTING_COEF_KEY";


MatrixBuilder::HashType GAJK::my_hash_(const string &,
                                     unsigned int deriv,
                                     const Wavefunction &wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto eri_key=options().get<string>(ERI_opt);
    auto eri=create_child_from_option<FourCenterIntegral>(ERI_opt);
    auto D_hash=bphash::hash_to_string(wfn.opdm->my_hash());
    return eri_key+D_hash+eri->my_hash(deriv,wfn,bs1,bs2,bs1,bs2);
    return "";
}

MatrixBuilder::ReturnType GAJK::calculate_(const string &,
                          unsigned int deriv,
                          const Wavefunction & wfn,
                          const BasisSet &bs1,
                          const BasisSet &)
{
    //Check if we were supposed to get the result from the cache
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache, but FORCE_CACHE=true");

    //Get the alpha density matrix
    auto sharedD=convert_to_GA(*wfn.opdm->get(Irrep::A,Spin::alpha));
    const auto& D=*sharedD;

    //Make and initialize a module that will provide the ERIs
    //Should really use bs2 and not assume symmetric...
    const auto Ints=create_child_from_option<FourCenterIntegral>(ERI_opt);
    Ints->initialize(deriv,wfn,bs1,bs1,bs1,bs1);

    //J, K, and D dimensions
    const size_t nbf=bs1.n_functions();
    const std::array<size_t,2> dims({nbf,nbf});

    //Will be J and K
    GATensor J(dims,0.0), K(dims,0.0);

    //A module that can form the norm of the "diagonal" elements of the ERIs
    auto schwarz_metric=create_child<MatrixBuilder>("PSR_Sieve");

    //These are said "diagonal" elements
    const auto& metric=
         *convert_to_eigen(*schwarz_metric->calculate("",deriv,wfn,bs1,bs1)[0]);

//    //This is the actual screen object
//    GATensorImpl<2> temp_D(D);
//    SchwarzScreen sieve(metric,temp_D,bs1);

    //This is an iterator over the unique shell quartets
    ShellQuartetItr quarts(bs1);

    //Pointers to the local data
    auto Jdata=J.my_data(),Kdata=K.my_data(),Ddata=D.my_data();
    auto pJ=Jdata.data();
    auto pK=Kdata.data();
    auto pD=Ddata.data();

    while(quarts)
    {
        const auto& shell=*quarts;//Current shell

//        if(!sieve.is_good(shell))//See if shell quartet gets sieved out
//        {
//            ++quarts;
//            continue;
//        }

        //Factor to scale the integral (see JK_build notes in dox/)
        const double total_deg=quarts.degeneracy();

        //Calculate the ERI
        const double* buffer=
                Ints->calculate(shell[0],shell[1],shell[2],shell[3]);

        //Digest it
        size_t counter=0;
        for(const auto& bf_quart:quarts)
        {
            const double value=buffer[counter++]*total_deg;
            const size_t mu=bf_quart[0],nu=bf_quart[1],
                         lambda=bf_quart[2],sigma=bf_quart[3];

            pJ[mu*nbf+nu]+=pD[lambda*nbf+sigma]*value;
            pJ[lambda*nbf+sigma]+=pD[mu*nbf+nu]*value;
            pK[mu*nbf+lambda]-=pD[nu*nbf+sigma]*value;
            pK[mu*nbf+sigma]-=pD[nu*nbf+lambda]*value;
            pK[nu*nbf+sigma]-=pD[mu*nbf+lambda]*value;
            pK[nu*nbf+lambda]-=pD[mu*nbf+sigma]*value;
        }
        ++quarts;
    }

    //Update GA
    J.set_value(J.my_shape(),pJ);
    K.set_value(K.my_shape(),pK);

    //Symmetrize
    GATensor J_final=symmetrize(J),K_final=symmetrize(K);

    auto rvJ=make_shared<GATensorImpl<2>>(J_final);
    auto rvK=make_shared<GATensorImpl<2>>(K_final);

    //Return following Pulsar's convention
    return {rvJ,rvK};
}

}
