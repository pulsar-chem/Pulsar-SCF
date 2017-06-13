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
    //A core guess is equivalent to guessing a density of 0 and running a single
    //SCF iteration.  Is there a clever way to exploit that in the SCF iteration
    //so we don't even have this function anymore?

    //Get basis set, occupations, core hamiltonian, and S
    //See comments in SCF.cpp for more information
    const auto bs=wfn.system->get_basis_set("PRIMARY");
    const auto temp_occs=guess_occ(wfn)->get(Irrep::A,Spin::alpha);
    const auto occs=*convert_to_eigen(*temp_occs);
    const size_t nocc=std::accumulate(occs.data(),occs.data()+occs.size(),0);

    auto H_builder=create_child_from_option<MatrixBuilder>("H_KEY");
    const auto H_temp=H_builder->calculate("",deriv,wfn,bs,bs)[0];
    const auto H=*convert_to_eigen(*H_temp);

    const auto S_builder=create_child_from_option<MatrixBuilder>("S_KEY");
    const auto S_temp=S_builder->calculate("",deriv,wfn,bs,bs)[0];
    const auto S=*convert_to_eigen(*S_temp);

    Eigen::GeneralizedSelfAdjointEigenSolver<matrix_type> eSC(H, S);
    matrix_type C = eSC.eigenvectors().leftCols(nocc);
    const matrix_type D = C*C.transpose();

    return {update_wfn(wfn,&C,&D,&occs),{0.0}};
}

}//End namespace
