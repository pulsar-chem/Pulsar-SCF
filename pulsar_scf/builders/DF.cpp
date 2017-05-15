#include "pulsar_scf/builders/DF.hpp"
#include <pulsar/modulebase/TwoCenterIntegral.hpp>
#include <pulsar/modulebase/ThreeCenterIntegral.hpp>
#include "pulsar_scf/MatrixFillFxns.hpp"

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using tensor_type=EigenTensorImpl<3>::tensor_type;
using tensor_matrix=Eigen::TensorMap<Eigen::Tensor<const double,2>>;

namespace pulsar_scf {

MatrixBuilder::ReturnType Metric::calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto MetricInts=
            create_child_from_option<TwoCenterIntegral>("METRIC_INTS_KEY");
    MetricInts->initialize(deriv,wfn,bs1,bs2);
    const bool is_symmetric=bs1==bs2;
    auto Metric(std::move(is_symmetric ?
                       detail_::fill_symmetric<matrix_type>(MetricInts,bs1) :
                       detail_::fill_asymmetric<matrix_type>(MetricInts,bs1,bs2)
           ));
    return {std::make_shared<EigenMatrixImpl>(std::move(Metric))};
}

Rank3Builder::ReturnType DFInts::calculate_(const std::string &,
                      unsigned int deriv,
                      const pulsar::Wavefunction & wfn,
                      const pulsar::BasisSet & bs1,
                      const pulsar::BasisSet & bs2,
                      const pulsar::BasisSet & bs3)
{
    auto DFInts=
            create_child_from_option<ThreeCenterIntegral>("DF_INTS_KEY");
    DFInts->initialize(deriv,wfn,bs1,bs2,bs3);
    const bool is_symmetric=bs1==bs2;
    const size_t nbf1=bs1.n_functions(),
                 nbf2=bs2.n_functions(),
                 nbf3=bs3.n_functions();
    const size_t nshell1=bs1.n_shell(),
                 nshell2=bs2.n_shell(),
                 nshell3=bs2.n_shell();
    Eigen::Tensor<double,3> ints(nbf1,nbf2,nbf3);
    for(size_t i=0;i<nshell1;++i)
    {
        const size_t nbfi=bs1.shell(i).n_functions();
        const size_t offi=bs1.shell_start(i);
        for(size_t j=0;j<nshell2;++j)
        {
            const size_t nbfj=bs2.shell(j).n_functions();
            const size_t offj=bs2.shell_start(j);
            for(size_t k=0;k<nshell3;++k)
            {
                const size_t nbfk=bs3.shell(k).n_functions();
                const size_t offk=bs3.shell_start(k);
                auto buffer=DFInts->calculate(i,j,k);
                if(buffer==nullptr)continue;
                for(size_t n1=0,counter=0;n1<nbfi;++n1)
                    for(size_t n2=0;n2<nbfj;++n2)
                        for(size_t n3=0;n3<nbfk;++n3,++counter)
                            ints(n1+offi,n2+offj,n3+offk)=buffer[counter];

            }
        }
    }

    return {std::make_shared<EigenTensorImpl<3>>(std::move(ints))};
}

Rank3Builder::ReturnType DFCoef::calculate_(const std::string&,
                                            unsigned int deriv,
                                            const Wavefunction& wfn,
                                            const pulsar::BasisSet& dfbs,
                                            const pulsar::BasisSet& bs1,
                                            const pulsar::BasisSet& bs2)
{
    const size_t ndf=dfbs.n_functions();
    const auto metric_ints=
            create_child_from_option<MatrixBuilder>("METRIC_KEY");
    const auto df_ints=
            create_child_from_option<Rank3Builder>("DF_INTS_KEY");
    auto Jmetric =
            *convert_to_eigen(*metric_ints->calculate("",deriv,wfn,dfbs,dfbs)[0]);
    matrix_type Linv_temp=
          Jmetric.llt().matrixL().solve(matrix_type::Identity(ndf,ndf));
    tensor_type d_Qls=
          *convert_to_eigen(*df_ints->calculate("",deriv,wfn,dfbs,bs1,bs2)[0]);
    tensor_matrix Linv(Linv_temp.data(), ndf,ndf);
    std::array<Eigen::IndexPair<int>,1> idx({Eigen::IndexPair<int>(1,0)});
    return {std::make_shared<EigenTensorImpl<3>>(std::move(Linv.contract(d_Qls,idx)))};
}

}//End namespace
