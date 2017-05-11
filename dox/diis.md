DIIS                                                                {#DIIS}
====

This page describes DIIS (direct inversion of the iterative subspace) as it is
implemented in this project.

## Theory

In an iterative procedures one recomputes some quantity repeatedly until
the quantity is converged.  Let us assume that on iteration \f$i\f$ this
quantity is \f$x_i\f$.  After \f$m\f$ iterations we have amassed \f$m\f$ such
"guesses".  Let us assume there exists some function \f$f\f$ such that for the
true answer, \f$x\f$, \f$f(x)=0\f$ (for \f$f(x)=c\ c\neq0\f$ this can be
trivially accomplished by defining \f$f^\prime (x)=f(x)-c\f$ and instead work
with \f$f^\prime(x)\f$).  This means for each \f$x_i\f$ we can compute an error
(note an error is the difference between a computed value and the true value, a
residual is the difference betweeen a computed value and a model value; 0 is the
true value...)  \f$e_i=f(x_i)-0=f(x_i)\f$.

Each of our guesses can be thought of as \f$x_i=x-\Delta x_i\f$ where
\f$\Delta x_i\f$ is the part of \f$x\f$ missing from \f$x_i\f$.  Without knowing
\f$x\f$ we do not know \f$\Delta x_i\f$.  Instead of using \f$x\f$ we can use a
model for \f$x\f$, \f$x^\prime\f$.  The result is a residual,
\f$\Delta x_i^\prime=x^\prime-x_i\f$. For our model we assume that
\f$x^\prime\f$ can be written as a linear combination of our previous guesses,
that is:
\f[
x^\prime =\sum_{i=1}^{m}c_i x_i
\f]
To determine how good our model is we can apply \f$f\f$.  The result is:
\f[
f(x^\prime)=\sum_{i=1}^{m}f(c_ix_i)=\sum_{i=1}^mc_if(x_i)=\sum_{i=1}^{m}c_ie_i
\f]
where we have assumed \f$f\f$ is linear in \f$x_i\f$.  Thus the error (not the
residual) in our fit is given by a linear combination of the errors in our
previous guesses.  Thus the best guess is obtained by choosing the \f$c\f$'s to
minimize \f$f(x^\prime)\f$.  In general, however, \f$f\$ may have both maxima
and minima so taking the derivative and setting it to zero does not necessarilly
lead to a minimum.  This can trivially be circumvented by instead minimizing the
square of \f$f(x^\prime)\f$:
\f[
f(x^\prime)^2=\sum_{i=1}^{m}\sum_{j=1}^mc_i^\dagger \langle e_i, e_j\rangle c_j
\f]
Assuming the coefficients are real, and defining a matrix \f$B\f$ such that
\f$B_{ij}=\langle e_i, e_j\rangle\f$ this leads to:
\f[
\frac{\partial f(x^\prime)^2}{\partial c_k}=
   \sum_{i=1}^{m}c_iB_{ij}+
   \sum_{j=1}^m B_{ij}c_j=
   2\sum_{i=1}^{m}c_iB_{ij}
\f]

There is one final issue to consider in the DIIS procedure, normalization.  If
we assume that each \f$x_i\approx x\f$, then:
\f[
x^\prime =\sum_{i=1}^m c_i x_i\approx x\sum_{i=1}^m c_i
\f]
That is we get back approximately \f$\sum_{i=1}^m\f$ times \f$x\f$.  Hence we
need to insist that \f$\sum_{i=1}^mc_i=1\f$.  This is most easily done via
Lagrange's method of undertermined multipliers.  In such a case we define a
function:
\f[
g(x^\prime)=f(x^\prime)^2+\lambda(1-\sum_{i=1}^m c_i),
\f]
which for all values of \f$x^\prime\f$ is equal to \f$f(x^\prime)^2\f$; however,
the minimum is now at:
\f[
\frac{\partial g(x^\prime)}{\partial c_k}=
\frac{\partial f(x^\prime)^2}{\partial c_k}-\lambda=
\sum_{i=1}^{m}c_iB_{ij}-\lambda
\f]
Where in the last line we have absorbed the 2 into lambda.

## General Considerations for Implementation

Many codes, manuals, etc. discuss when to turn on DIIS.  This is, in our
opinion, a moot point from the perspective of this class.  If you don't want to
do DIIS yet don't call this class...In a similar vein, DIIS requires at least
two guesses in order to extrapolate (if you only have one, the only possible
weight for that guess is 1, thus returning your guess).  Thus update is only
called starting with the second guess.

Another interesting question is whether to use the resulting guess to check for
convergence.  More specifically on each iteration of an iterative procedure, one
generates a guess \f$x_i\f$ and its associated error \f$e_i\f$.  DIIS then
extrapolates this guess to \f$x\f$, which has associated error \f$e\f$.  Now
should the iteration use \f$e_i\f$ or \f$e\f$ to check for convergence?

## DIIS class

DIIS is handled via the DIIS class.  This class is designed to work with any
form of data structure via a visitor pattern design.  The first template type
parameter, `T`, is the type of your data structure.  The second template type
parameter is the type of a functor to be used for the dot prouduct of two
instances of type `T`.  Assuming `T=Eigen::MatrixXd` we have provided an
implementation for you:
~~~{.cpp}
struct eigen_dot
{
    double operator()(const Eigen::MatrixXd& v1, const Eigen::MatrixXd& v2)const
    {

        return v1.cwiseProduct(v2).sum();
    }
};
~~~
Otherwise you will need to provide a functor type yourself (replace
`Eigen::MatrixXd` with your `T`).  Similarly, the DIIS procedure also needs
a way of taking the linear combination.  A functor capable of doing this is
provided by third template type parameter.  Again, if you use `Eigen::MatrixXd`
an implementation is provided:
~~~{.cpp}
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
~~~
If you are not using `T=Eigen::MatrixXd` then again you will need to provide
your own functor (again replace `Eigen::MatrixXd` with the type of your data
structure).  Users should note that the Lagrange multiplier is the last
coefficient in the list passed to this functor.
