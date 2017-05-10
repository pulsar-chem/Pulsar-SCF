#pragma once
#include<pulsar/math/EigenImpl.hpp>
#include<deque>

namespace pulsar_scf {
namespace detail_ {
///Functor for dotting two Eigen matrices
struct eigen_dot
{
    double operator()(const Eigen::MatrixXd& v1, const Eigen::MatrixXd& v2)const
    {

        return v1.cwiseProduct(v2).sum();
    }
};

///Functor for taking linear combinations of Eigen matrices
struct eigen_combiner
{
    Eigen::MatrixXd operator()(const double* cs,
                               const std::deque<Eigen::MatrixXd>& vs)const
    {
        Eigen::MatrixXd result=cs[0]*vs[0];
        for(size_t i=1;i<vs.size();++i)
            result+=cs[i]*vs[i];
        return result;
    }
};
}

/** \brief A helper class for performing DIIS extrapolation
 *
 *  For theoretical details pertaining to DIIS see [here](#diis).
 *
 *  \note  Although queue is more natural than deque the
 *  stupidity of the STL in not letting us use the index operator for queue is
 *  why we use deque...
 *
 *  The class assumes that B is relatively small (and hence also the resulting
 *  linear expansion coefficients) in order to hide the user from the fact
 *  that internally we use Eigen.
 *
 * \tparam T the type of the matrix
 */
template<typename T,
         typename dot_product=detail_::eigen_dot,
         typename axpy_t=detail_::eigen_combiner>
class DIIS{
private:
    dot_product dp_;///<Functor for performing the dot product
    axpy_t axpy_;///<Functor for linear combination of the previous guesses
    std::deque<T> vectors_;///<Previously found guesses
    std::deque<T> errors_;///<Previously found errors
    const size_t max_;///<Maximum number of previous guesses to maintain
    Eigen::MatrixXd B_;///<Inner product of guesses

    ///Updates the last row and column of B
    void update_B(const T& ei){
        //Index of current error noting errors already includes ei
        const size_t i=errors_.size()-1;
        for(size_t j=0;j<i;++j)
        {
          const double value=dp_(ei,errors_[j]);
          B_(i,j)=value;
          B_(j,i)=value;
        }
        B_(i,i)=dp_(ei,ei);
    }

    ///Returns the new expansion coefficients using the current B matrix
    Eigen::VectorXd solve()const
    {
        const size_t dim=errors_.size();

        Eigen::MatrixXd A=Eigen::MatrixXd::Zero(dim+1,dim+1);
        A.topLeftCorner(dim,dim)=B_.topLeftCorner(dim,dim);
        A.row(dim).setConstant(-1.0);
        A.col(dim).setConstant(-1.0);
        A(dim,dim)=0.0;
        Eigen::VectorXd b=Eigen::VectorXd::Zero(dim+1);
        b(dim)=-1.0;
        return A.colPivHouseholderQr().solve(b);
    }

    ///Removes the oldest guess/error to stay within the requested number of
    void prune()
    {
        vectors_.pop_front();
        errors_.pop_front();
        Eigen::MatrixXd tempB=B_.bottomRightCorner(max_-1,max_-1);
        tempB.conservativeResize(max_,max_);
        B_=tempB;
    }


public:
    ///Given an initial guess and err makes a class capable of performing DIIS
    DIIS(const T& vi,const T& ei,size_t max_vec=5):
        vectors_({vi}),errors_({ei}),
        max_(max_vec),B_(Eigen::MatrixXd::Zero(max_vec,max_vec))
    {
        B_(0,0)=dp_(ei,ei);
    }

    /** \brief Updates the state of the class and returns the new guess
     *
     */
    T update(const T& vi, const T& ei){
        vectors_.push_back(vi);
        errors_.push_back(ei);
        if(vectors_.size()==max_+1)prune();
        update_B(ei);
        auto cs=solve();
        return axpy_(cs.data(),vectors_);
    }
};


}
