#include "pulsar_scf/builders/Builders.hpp"
#include "pulsar_scf/MatrixFillFxns.hpp"
#include <memory>

using namespace pulsar;
using namespace std;
using matrix_type=EigenMatrixImpl::matrix_type;
using ReturnType=MatrixBuilder::ReturnType;

namespace pulsar_scf {

//The key for the overlap matrix builder
const string S_opt="S_KEY";

//Implements the hash for an orthogonalizer
MatrixBuilder::HashType Orthogonalizer::my_hash_(const string &,
                                                  unsigned int deriv,
                                                  const Wavefunction &wfn,
                                                  const BasisSet &bs1,
                                                  const BasisSet &bs2)
{
    //An orthogonalizer only depends on the underlying overlap matrix builder
    auto S_key=options().get<string>(S_opt);
    auto Sbuilder=create_child_from_option<MatrixBuilder>(S_opt);
    return S_key+Sbuilder->my_hash("",deriv,wfn,bs1,bs2);
}


//This is flat out stolen from libint2's hartree-fock++.cc routine
ReturnType Orthogonalizer:: calculate_(const std::string &,
                                     unsigned int deriv,
                                     const Wavefunction & wfn,
                                     const BasisSet &bs1,
                                     const BasisSet &bs2)
{
    //Hack to ensure that cacheing occurs
    if(options().get<bool>("FORCE_CACHE"))
        throw PulsarException("Value not in cache and FORCE_CACHE = true");

    //TODO: the threshold should be made into an option
    const double threshold = 1e-6;
    //TODO: the abilitiy to do symmetric orthogonalization should become another
    //      module
    const bool symmetric=false;

    //Get an S builder, and then get S
    auto S_builder=create_child_from_option<MatrixBuilder>(S_opt);
    const auto S_temp=S_builder->calculate("",deriv,wfn,bs1,bs2)[0];
    const auto S=*convert_to_eigen(*S_temp);

    //Need to form S^-1/2
    //Do this by diagonalization and then raising each element to -1/2 power

    //Get eigen values/vectors
    Eigen::SelfAdjointEigenSolver<matrix_type> eig_solver(S);
    auto s = eig_solver.eigenvalues();
    auto U = eig_solver.eigenvectors();

    //In order to avoid numeric instabilities we skip eigenvalues that are close
    //to zero as defined by threshold
    //Note: This algorithm exploits that the eigenvalues are in increasing order

    //This is the number of good eigenvalues
    size_t n_cond = 0;
    //This is the number of eigenvalues
    const size_t n=s.rows();
    //TODO: If the number of eigenvalues to throw away is usually small this
    //      would be better implemented by starting from 0 and finding the first
    //      value >=threshold (note that this is std::lower_bound).  Using an
    //      std algorithm moves us to iterators so one would have to ensure the
    //      following can be done with iterators (note there are reverse
    //      iterators for looping over a container backwards)

    //We now loop backwards over the container
    for(size_t i=n;i>0;--i)
    {
      if (s(i-1) >= threshold)//If this element is above threshold count it
          ++n_cond;
      else //Otherwise we're done
          break;
    }

    //Grab the good eigenvalues, take the square root, and set them down the
    //diagonal
    //(element-wise operations in Eigen can only be done on arrays)
    auto sigma = s.bottomRows(n_cond).array().sqrt();
    auto sigma_sqrt = sigma.matrix().asDiagonal();
    //Invert the matrix
    auto sigma_invsqrt = sigma.inverse().matrix().asDiagonal();

    //Grab only the eigenvectors associated with the good eigenvalues
    auto U_cond = U.block(0, n-n_cond, n, n_cond);
    //Form X and X inverse
    matrix_type X = U_cond * sigma_invsqrt;
    matrix_type Xinv = U_cond * sigma_sqrt;

    // These are the additional steps required by symetric orthogonalization
    if (symmetric) {
      X = X * U_cond.transpose();
      Xinv = Xinv * U_cond.transpose();
    }

    return {make_shared<EigenMatrixImpl>(move(X)),
            make_shared<EigenMatrixImpl>(move(Xinv))};
}

}//End namespace
