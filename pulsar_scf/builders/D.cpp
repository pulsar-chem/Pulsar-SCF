#include "pulsar_scf/builders/D.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include <memory>
#include <Eigen/Dense>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
namespace pulsar_scf {

DerivReturnType CoreDensity::deriv_(size_t deriv,
                                     const Wavefunction & wfn)
{
    const auto bs=wfn.system->get_basis_set("PRIMARY");
    auto occs=guess_occ(wfn);
    const auto H=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("H_KEY")
                    ->calculate("",deriv,wfn,bs,bs)[0]);
    const auto S=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("S_KEY")
                    ->calculate("",deriv,wfn,bs,bs)[0]);
    Eigen::GeneralizedSelfAdjointEigenSolver<matrix_type> eSC(H, S);
    matrix_type C = eSC.eigenvectors().leftCols(5);
    const matrix_type D = C*C.transpose();
    return {update_wfn(wfn,&C,&D),{0.0}};
}

}//End namespace
