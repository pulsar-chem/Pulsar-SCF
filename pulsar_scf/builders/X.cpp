#include "pulsar_scf/builders/X.hpp"
#include "pulsar_scf/MatrixFillFxns.hpp"
#include <pulsar/modulebase/OneElectronIntegral.hpp>
#include <memory>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;
namespace pulsar_scf {

///This is flat out stolen from libint2's hartree-fock++.cc routine

ReturnType Orthogonalizer:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    auto S=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("S_KEY")
                             ->calculate("",deriv,wfn,bs1,bs2)[0]);
    Eigen::SelfAdjointEigenSolver<matrix_type> eig_solver(S);
    auto U = eig_solver.eigenvectors();
    auto s = eig_solver.eigenvalues();
    const double threshold = 1e-6;
    const bool symmetric=false;
    size_t n_cond = 0;
    const size_t n=s.rows();
    for(size_t i=n;i>0;--i)
    {
      if (s(i-1) >= threshold)++n_cond;
      else break;
    }

    auto sigma = s.bottomRows(n_cond).array().sqrt();
    auto sigma_sqrt = sigma.matrix().asDiagonal();
    auto sigma_invsqrt = sigma.inverse().matrix().asDiagonal();

    // make canonical X/Xinv
    auto U_cond = U.block(0, n-n_cond, n, n_cond);
    matrix_type X = U_cond * sigma_invsqrt;
    matrix_type Xinv = U_cond * sigma_sqrt;
    // convert to symmetric, if needed
    if (symmetric) {
      X = X * U_cond.transpose();
      Xinv = Xinv * U_cond.transpose();
    }
    return {std::make_shared<EigenMatrixImpl>(std::move(X)),
            std::make_shared<EigenMatrixImpl>(std::move(Xinv))};
}

}//End namespace
