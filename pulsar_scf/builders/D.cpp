#include "pulsar_scf/builders/D.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include <memory>
#include <numeric>
#include <Eigen/Dense>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
namespace pulsar_scf {

DerivReturnType CoreDensity::deriv_(size_t deriv,
                                     const Wavefunction & wfn)
{
    const auto bs=wfn.system->get_basis_set("PRIMARY");
    const auto occs=*convert_to_eigen(*guess_occ(wfn)->get(Irrep::A,Spin::alpha));
    const size_t nocc=std::accumulate(occs.data(),occs.data()+occs.size(),0);
    const auto H=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("H_KEY")
                    ->calculate("",deriv,wfn,bs,bs)[0]);
    const auto S=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("S_KEY")
                    ->calculate("",deriv,wfn,bs,bs)[0]);
    Eigen::GeneralizedSelfAdjointEigenSolver<matrix_type> eSC(H, S);
    matrix_type C = eSC.eigenvectors().leftCols(nocc);
    const matrix_type D = C*C.transpose();
    return {update_wfn(wfn,&C,&D,&occs),{0.0}};
}

}//End namespace
