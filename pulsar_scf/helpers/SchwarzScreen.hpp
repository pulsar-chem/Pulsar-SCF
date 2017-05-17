#pragma once
#include<pulsar/modulebase/MatrixBuilder.hpp>

/** \brief
 *
 *
 *  This class works based off the Cauchy-Schwarz inequality:
 *  \f[
 *     ||\langle x,y \rangle||^2 \le \langle x,x\rangle \times \langle y,y\rangle
 *  \f]
 *  Basically, think of \f$x\f$ or \f$y\f$ as a shell pair, then if we pre-tabulate
 *  the value of \f$\langle x,x\rangle\f$ (which will be the same set as
 *  \f$\langle y,y\rangle\f$ owing to \f$x\f$ and \f$y\f$ being dummy indices)
 *  then we can estimate the magnitude of any mixed inner product by multiplying
 *  their values together.
 */
namespace pulsar_scf {

class SchwarzScreen: public pulsar::MatrixBuilder {
public:
    SchwarzScreen(ID_t id):MatrixBuilder(id){}
    ReturnType calculate_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
    HashType my_hash_(const std::string & key,
                          unsigned int deriv,
                          const pulsar::Wavefunction & wfn,
                          const pulsar::BasisSet & bs1,
                          const pulsar::BasisSet & bs2);
};

}
