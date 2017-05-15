# JK Builders                                         {#jkbuilders}
This page is meant to describe the various JK builders available and explain the
implementations.

## Theory

Particularly when describing the Hartree-Fock method there are two matrices that
take a central role, namely the Coloumb matrix, \f$J\f$, and the exchange
matrix, \f$K\f$.  The element \f$J_{\mu,\nu}\f$ describes the electron repulsion
of an electron in orbital \f$\mu\f$ with that of an electron in orbital
\f$\nu\f$.  The element \f$K_{\mu,\nu}\f$ describes the change in energy if
electron \f$\mu\f$ and \f$\nu\f$ simaltaneously tunnel.  Formulaically, in the
atomic orbital (AO) basis set, we have:
\f[
J_{\mu\nu}=\sum_{\lambda}^N\sum_{\sigma}^{N}D_{\lambda\sigma}(\mu\nu|\lambda\sigma)
\f]
and
\f[
K_{\mu\nu}=\sum_{\lambda}^N\sum_{\sigma}^ND_{\lambda\sigma}(\mu\sigma|\lambda\nu)
\f]
where \f$D\f$ is the one-particle density matrix in the \f$N\f$ orbital, AO
basis set and \f$(\mu\nu|\lambda\sigma)\f$ are the two-electron, four-center
electron repulsion integrals in chemist's notation *i.e.* the bra is the indices
for electron 1 and the ket is the indices for electron 2.  Note that for
unrestricted SCF \f$D\f$ in \f$J\f$ is \f$D^\alpha+D^\beta\f$, whereas that in
\f$K\F$ is \f$D^\alpha\f$ for the \f$\alpha\f$ block of the Fock matrix and
\f$D^\beta\f$ for the \f$\beta\f$ block.

The above are the fundamental equations for \f$J\f$ and \f$K\f$, but
programmatically they are rarely implemented as shown tensor contractions.  How
they are implemented tends to depend on how the integrals are formed and
comprises the next several sections. It is helpful to first review the symmetry
of the ERIs.  Specifically, note:
\f[
(\mu\nu|\lambda\sigma)=(\nu\mu|\lambda\sigma)=
(\mu\nu|\sigma\lambda)=(\nu\mu|\sigma\lambda)=
(\lambda\sigma|\mu\nu)=(\sigma\lambda|\mu\nu)=
(\lambda\sigma|\nu\mu)=(\sigma\lambda|\nu\mu).
\f]
These symmetries follow from the interchange of which electron is 1 and which is
2 (1 and 2 are just dummy indicies) as well as swapping the orbitals for a given
electron between the bra and ket (a symmetry only present if the orbitals are
real).

## Direct

Arguably the simplest integral implementation is what is called direct.  In a
direct scheme each integral (actually a block of integrals) is computed when it
is needed.  Typically these blocks are so called shell quartets (a shell is a
set of basis functions that have the same exponents, center, and contraction
coefficients).  We thus will have four nested loops over shells (one loop for
each index of the ERI).  We then additionally have four more nested loops to
digest the ERI (these four loops run over the basis functions in each shell).
In order to take advantage of the symmetries above for the integral
\f$(\mu\nu|\lambda\sigma)\f$ we first need to ensure that we are only computing
unique pairs \f$\mu\nu\f$ in the bra and unique pairs \f$\lambda\sigma\f$
in the ket (we do this by requiring \f$\mu\ge\nu\f$ and \f$\lambda\ge\sigma\f$).

This only utilizes the second, third, and fourth symmetries.  To use the
remainder we additionally have to ensure \f$\mu\nu\ge\lambda\sigma\f$. In
turn we have the following loop ranges:
\f[
s_1\in[0,nshells); s_2\in[0,s_1]; s_3\in[0,s_1]
\f]
The fourth loop is more complicated and is:
\f[
s_4\in[0,s_3]
\f]
if \f$s_1\neq s_3\f$, otherwise:
\f[
s_4\in[0,s_2]
\f]

As written above \f$J\f$ and \f$K\f$ are in terms of different ERIs.  We need to
write them in terms of the same ERI (or else we will have to recompute the shell
quartets for \f$J\f$ and \f$K\f$).  Keeping \f$J\f$ fixed and renaming dummy
indicies in \f$K\f$ so that the ERI indicies match gives (assuming the density
matrix is symmetric):
\f[
K_{\mu\sigma}=\sum_{\nu}^{N}\sum_{\lambda}^ND_{\nu\lambda}(\mu\nu|\lambda\sigma).
\f]
Next, we have to figure out which elements of \f$J\f$ and \f$K\f$ to update for
each integral by applying the symmetry relations.  Since the summation in
\f$J\f$ is over all permutations of \f$\lambda\sigma\f$ and not just pairs
less than \f$\mu\nu\f$ we miss updating:
\f[
J_{\lambda\sigma}=D_{\mu\nu}(\mu\nu|\lambda\sigma).
\f]
For \f$K\f$ we miss a similar term:
\f[
K_{\nu\lambda}=D_{\mu\sigma}(\mu\nu|\lambda\sigma),
\f]
but also miss cases in which \f$\nu>\mu\f$ and thus have two additional terms:
\f[
K_{\nu\sigma}=D_{\mu\lambda}(\mu\nu|\lambda\sigma)
\f]
and
\f[
K_{\mu\lambda}=D_{\nu\sigma}(\mu\nu|\lambda\sigma).
\f]

Finally, we have to account for some undercounting.  Note again that the
summations in \f$J\f$ and \f$K\f$ are not over pairs, but permutations.  Before
working out the undercountings, note that making matters more complicated, in
the inner four loops one typically runs all indices to completion, instead of
ensuring only unique pairs are considered.  If you are curious, this is a
vectorization optimization stemming from the fact that most integral libraries
always return the full quartet and that we thus can only access this integral
buffer contigiously if we visit all indices.  Thus, unless \f$\s_1=s_2\f$ we
will not include the pair \f$\nu\mu\f$.  Similarly, we
will not include the pair \f$\sigma\lambda\f$ unless \f$s_3=s_4\f$.  Finally,
if \f$s_1\neq s_3\f$ we miss the pair \f$\lambda\mu\f$ and unless
\f$s_1\eq s_3\f$ and \f$s_2\eq s_4\f$ we miss the pair \f$\sigma\nu\f$  Putting
that all together, for each of the following that are not true, multiply by 2:
\f$s_1\eq s_2\f$, \f$s_3\eq s_4\f$, and \f$s_1\eq s_3 \wedge s_2\eq s_4\f$.

## Density Fit

Another popular way to build J and K is with density fitting.  The idea is that,
at least for Coulomb four-center, two-electron integral the integral can be
interpreted as the interaction of electron density assoicated with two atomic
orbitals.  If we expand each density in terms of an \f$M\f$ dimensional
auxillary basis set \f$|Q)\f$ the result will contain less indices.  Given a
density \f$\rho_{\mu\nu}\f$ this leads to
\f[
\rho_{\mu\nu}\approx\sum_{Q}^Md^Q_{\mu\nu}|Q).
\f]
In general the dimension of \f$|Q)\f$ is much greater than that of the atomic
orbital basis set and we have an overdetermined system of equations.  One can
then use any fitting procedure one likes, but traditionally the coefficients are
choosen such that:
\f[
d^Q_{\mu\nu}=\sum_{P}^M(\mu\nu|P)\left[\mathbf{J}^{-1}\right]_{PQ}
\f]
It can be shown (TODO show it) that this choice minimizes the error in the
electric field generated by the approximate \f$\rho_{\mu\nu}\f$.


In this new basis set, the integrals are now:
\f{eqnarray*}{
(\mu\nu|\lambda\sigma)&=&\sum_{P}^M\sum_{Q}^M(\mu\nu|P)[\mathbf{J}^{-1}]_{PQ}(Q|\lambda\sigma)\nonumber\\
&=&\sum_{P}^M(\mu\nu|P)T^{P}_{\lambda\sigma}\\
T^{P}_{\lambda\sigma}&\equiv&\sum_{Q}^M[\mathbf{J}^{-1}]_{PQ}(Q|\lambda\sigma)
\f}
and our \f$J\f$ and \f$K\f$ builds are now:
\f[
J_{\mu\nu}=\sum_{\lambda}^N\sum_{\sigma}^{N}\sum_{P}^MD_{\lambda\sigma}(\mu\nu|P)T^{P}_{\lambda\sigma}
\f]
and
\f[
K_{\mu\nu}=\sum_{\lambda}^N\sum_{\sigma}^N\sum_{P}^MD_{\lambda\sigma}(\mu\sigma|P)T^{P}_{\lambda\nu}
\f]

The traditional \f$J\f$ and \f$K\f$ builds cost \f$\mathcal{O}(N^4)\f$ as written
the density-fit ones actually cost more \f$\mathcal{O}(N^4M)\f$.  However, note
that if we use the molecular orbital coefficients instead of the density matrix,
\f$J\f$ and \f$K\f$ can be written as:
\f{eqnarray*}{
J_{\mu\nu}&=&\sum_{\lambda}^N\sum_{\sigma}^{N}\sum_{P}^M\sum_i^{n}
            C_{\lambda i}C^*_{i\sigma}(\mu\nu|P)T^{P}_{\lambda\sigma}\nonumber\\
          &=&\sum_{\sigma}^{N}\sum_{P}^M\sum_i^{n}C^*_{i\sigma}(\mu\nu|P)T^{P}_{i\sigma}
\f}
which scales as \f$\mathcal{O}(N^3Mn)\f$ and for \f$K\f$:
\f{eqnarray*}{
K_{\mu\nu}&=&\sum_{\lambda}^N\sum_{\sigma}^N\sum_{P}^M\sum_i^{n}
            C_{\lambda i}C^*_{i\sigma}(\mu\sigma|P)T^{P}_{\lambda\nu}\nonumber\\
          &=&\sum_{P}^M\sum_i^{n}T^{*P}_{\mu i}T^{P}_{i\nu}
\f}
which scales as \f$\mathcal{O}(N^2Mn)\f$.  The formation of the common
intermediate\f$T^{P}_{i\nu}\f$ costs \f$\mathcal{O}(N^2Mn)\f$ which makes
forming $K$ cheaper than in direct.  If we break \f$J\f$ into two pieces:
\f[
T^P\equiv\sum_{\sigma}^N\sum_i^nC^*_{i\sigma}T^{P}_{i\sigma}
\f]
and
\f[
J_{\mu\nu}=(\mu\nu|P)T^P
\f}
the formation of \f$T^{P}_{i\sigma}\f$ dominates the building of \f$J\f$ and we
obtain an overall scaling reduction with respect to direct.

## Core

TODO: expand me

Basically use same trick as density fit to get one index to occupied.  Unlike
DF we have to store a rank 4 tensor, but if we can do it we get a savings and
our final cost is \f$\mathcal{O}(N^3n)\f$ (recall \f$N<M\f$).

## Disk

TODO: expand me
