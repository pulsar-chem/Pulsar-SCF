#pragma once
#include <pulsar/datastore/Wavefunction.hpp>
#include <pulsar/math/IrrepSpinMatrix.hpp>
#include <pulsar/modulebase/All.hpp>
namespace pulsar_scf {

///This function returns the number of shells in a Pulsar basis set accounting
///for the fact that some of them may be general shells
size_t nshells(const pulsar::BasisSet& bs);


/** \brief This function guesses occupations for our SCF computation
 *
 *  If `wf.occupations` is set they are simply returned.  Otherwise, if the
 *  number of electrons is even and the multiplicity is singlet all orbitals are
 *  assumed to be doubly occupied.  I haven't worried about other guesses yet...
 */
std::shared_ptr<const pulsar::IrrepSpinVectorD> guess_occ(const pulsar::Wavefunction& wf);

pulsar::Wavefunction update_wfn(const pulsar::Wavefunction& wf,
                                const Eigen::MatrixXd *Ca,
                                const Eigen::MatrixXd *Da,
                                const Eigen::VectorXd *Epsa,
                                const Eigen::MatrixXd *Cb=nullptr,
                                const Eigen::MatrixXd *Db=nullptr,
                                const Eigen::VectorXd *Epsb=nullptr
                                );

///Code factorization for MatrixBuilder implementations
pulsar::MatrixBuilder::ReturnType
           matrix_builder_kernel(unsigned int deriv,
                                 const pulsar::Wavefunction& wf,
                                 const pulsar::BasisSet& bs1,
                                 const pulsar::BasisSet& bs2,
                                 pulsar::TwoCenterIntegral& Ints);

}//End namespace
