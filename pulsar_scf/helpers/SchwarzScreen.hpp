#pragma once
#include<pulsar/modulebase/MatrixBuilder.hpp>
#include<pulsar/math/TensorImpl.hpp>
#include<pulsar/math/MinMax.hpp>
#include "pulsar_scf/PulsarSCF.hpp"

namespace pulsar_scf {

/**
*  The results of this builder are a matrix \f$S_{ij}=\sqrt{(ij|ij)^2}\f$.
*/
MATRIX_BUILDER(SchwarzMetric)


/** \brief This is the actual object to use to screen integrals
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
 *
 *  At the end of the day a small integral can still have a large impact on the
 *  resulting J or K matrix depending on the density matrix, thus the final
 *  decision is made by multiplying the left and right of the Cauchy-Schwarz
 *  inequality by the norm of the shell block of the density
 */
class SchwarzScreen{
   const Eigen::MatrixXd& metric_;
   const double threshold_;
   Eigen::MatrixXd density_;
public:
   SchwarzScreen(const Eigen::MatrixXd& metric_,
                 const Eigen::MatrixXd& density,
                 const pulsar::BasisSet& bs,
                 double threshold_=std::numeric_limits<double>::epsilon());
   //Converts from generic density matrices
   SchwarzScreen(const Eigen::MatrixXd& metric_,
                 const pulsar::TensorImpl<2,double>& density,
                 const pulsar::BasisSet& bs,
                 double threshold_=std::numeric_limits<double>::epsilon()):
       SchwarzScreen(metric_,*pulsar::convert_to_eigen(density),bs,threshold_)
   {}

   ///Returns true if this quartet is greater than the threshold
   bool is_good(const std::array<size_t,4>& quartet)const
   {
       //These are the blocks touched by this quartet
       const double max_d= pulsar::max(
                   density_(quartet[0],quartet[1]),
                   density_(quartet[2],quartet[3]),
                   density_(quartet[0],quartet[2]),
                   density_(quartet[0],quartet[3]),
                   density_(quartet[1],quartet[2]),
                   density_(quartet[1],quartet[3])
       );
       return max_d*metric_(quartet[0],quartet[1])*
               metric_(quartet[2],quartet[3])>threshold_;
   }
};


}
