#include "pulsar_scf/SCF.hpp"
#include "pulsar_scf/helpers/DIIS.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include <numeric>

using namespace pulsar;

//Typedef to save me from writing this out everytime
using matrix_type=EigenMatrixImpl::matrix_type;

//Computes the nuclear/nuclear repulsion of the system contained in the wf
//TODO: Make this into an EnergyMethod (lots of ways to do this better, want to
//      be able to take advantage of them eventually)
double nuclear_nuclear_repulsion(const Wavefunction& wfn)
{
    double ENuc=0.0;
    //Note atoms reside in wfn.system as a pointer
    for(auto ai: *wfn.system)
    {
        for(auto aj: *wfn.system)
        {
            if(aj==ai)break;//Prevents aj>=ai, i.e. we loop over pairs j<i
            const double dx=ai[0]-aj[0],
                         dy=ai[1]-aj[1],
                         dz=ai[2]-aj[2];
            ENuc+=ai.Z*aj.Z/std::sqrt(dx*dx+dy*dy+dz*dz);
        }
    }
    return ENuc;
}

namespace pulsar_scf {

//TODO: It will likely be somewhat challanging to do until one is very
//      comfortable with the current code, but UHF can actually be implemented
//      with only a few additional lines of code.  Hint you'll need to modify
//      the JK builder, the G builder, the F builder, and this file.  For the
//      G and F builder you'll only need to ensure there's a loop over spins
//      (and possibly tweak the scalars in G).  For this file the energy
//      equation needs tweaked.  JK is the biggest changes.
//
//      With UHF implemented I think it would be pretty easy to implement ROHF,
//      but my understanding of ROHF is a lot more fuzzy...


//TODO: This is a long term TODO.  Very fine control could be achieved by making
//      one entire SCF iteration a function call.  Under this scheme one could
//      switch out the F updater and C maker on a per iteration basis.
//      Furthermore, it could be done without messy flow control logic inside
//      the SCF iteration itself.  If one does this I foresee making a single
//      SCF iteration an EnergyMethod (after all each iteration generates an
//      energy) and having the loop itself be an EnergyMethod above that.
//      This also allows for more customization on the user's part for
//      determining convergence.


DerivReturnType SCF::deriv_(size_t deriv,
                            const Wavefunction & wfn)
{

    //Basis sets are put on system w/ key,value system
    //By convention "PRIMARY" is the key for the "usual" basis set
    //(i.e. the non-fitting one)
    auto bs=wfn.system->get_basis_set("PRIMARY");
    const double e_conv=1e-6, comm_conv=1e-6;
    const size_t max_iters=100;

    //TODO: which key to use for the basis, the definition of a converged SCF
    //      energy, the defintion of a converged commutator, and the maximum
    //      number of iterations should all be made into options.


    /* Throughout this section, where we are getting tensors, we need to be
     * careful as the returns are shared_ptr's to the tensors.  If we don't save
     * the returns to a temporary, that shared_ptr will go out of scope and our
     * references may be invalidated.  This is particularly relevant after the
     * calls to convert_to_eigen because if the argument to that function is not
     * actually an Eigen matrix, the downcast will fail, and new memory will be
     * allocated.
     */


    //Tensors in the wavefunction are indexed by irrep and spin; here we're
    //assuming C1 symmetry and RHF so we get the occupations for irrep A and
    //alpha spin
    auto temp_occs=guess_occ(wfn)->get(Irrep::A,Spin::alpha);

    //Tensors are stored in a type-erased format (basically as void*) we thus
    //need to cast them to whatever type we want to work with, here that's Eigen
    //Lines like convert_to_eigen will go away in a future version once either a
    //tensor library is settled upon or my TensorWrapper library is further
    //along
    const auto& occs=*convert_to_eigen(*temp_occs);

    //We assume that each occupation is either 0 or 1 and sum them up to get the
    //number of occupied orbitals
    const size_t nocc=std::accumulate(occs.data(),occs.data()+occs.size(),0);

    const double ENuc=nuclear_nuclear_repulsion(wfn);
    double oldE=ENuc;

    //Declare a pointer t a DIIS object, use a pointer so that we can use a
    //nullptr value to indicate that DIIS hasn't started
    std::unique_ptr<DIIS<matrix_type>> diis;

    //This gives us a function that is capable of building the Core Hamiltonian
    //This function is polymorphic to allow for H to be built in several ways
    auto H_builder=create_child_from_option<MatrixBuilder>("H_KEY");

    //Now we call that function, getting back a type erased tensor (the zero is
    //an artificat of the API in that we actually get back an array of tensors)
    const auto H_temp=H_builder->calculate("",deriv,wfn,bs,bs)[0];
    const auto H=*convert_to_eigen(*H_temp);


    //These lines are identical to H, except now for the overlap matrix
    auto S_builder=create_child_from_option<MatrixBuilder>("S_KEY");
    const auto S_temp=S_builder->calculate("",deriv,wfn,bs,bs)[0];
    const auto S=*convert_to_eigen(*S_temp);

    //Next we need the orthogonalizer, X, unlike H and S we actually get back
    //two orthogonalizers X and X^-1 (respectively in that order)
    const auto Xs=create_child_from_option<MatrixBuilder>("X_KEY")->
                       calculate("",deriv,wfn,bs,bs);
    const auto X=*convert_to_eigen(*Xs[0]);

    //TODO: This X is canonical orthoganilization.  Symmetric orthogonalization
    //      can be implemented somewhat trivially by making a MatrixBuilder that
    //      calls the canonical orthogonalizer and then right multiplies by
    //      the eigen vectors of S.  As X.cpp is setup, this could be done via
    //      an option, but this introduces additional flow logic.  This "if"
    //      statement can be avoided by polymorphism (well technically the if
    //      statement just occurs at the driver level...)


    //Finally we get the C1 symmetry alpha density from the wavefunction
    matrix_type D=*convert_to_eigen(*wfn.opdm->get(Irrep::A,Spin::alpha));

    //Copy the wavefunction to get one that is suitably allocated to be used as
    //a return
    Wavefunction newwfn(wfn);


    size_t iter=0;

    //Similar to H, S, and X there are multiple ways to compute the Fock
    //matrix so we get the function that we presently will use
    auto F_builder=create_child_from_option<MatrixBuilder>("F_KEY");

    //The actual SCF loop
    do{

        //Get our current F matrix
        const auto F_temp=F_builder->calculate("",deriv,newwfn,bs,bs)[0];
        auto F=*convert_to_eigen(*F_temp);

        //Compute the energy of our new F
        const double newE=D.cwiseProduct(H+F).sum()+ENuc;

        //This is the change in energy
        const double dE=newE-oldE;

        //This is the two-norm of the commutator
        const auto comm=(F*D*S-S*D*F);
        const double rmsComm=comm.norm();


        //TODO: Turn updating the Fock matrix into a MatrixBuilder.  This
        //      allows for polymorphic control over the update step, i.e.
        //      switching from say DIIS to second order
        if(!diis)//On first iteration just make our DIIS instance
            diis=std::make_unique<DIIS<matrix_type>>(F,comm);
        else //Otherwise let DIIS extrapolate F
            F=diis->update(F,comm);

        //TODO: Turn the generation of the new C matrices into a MatrixBuilder
        //      so it
        //      can be done polymorphically, say via NOGA or purification.
        //      In my head, this is done by combining the Fock builder and the
        //      Fock updater into one module that given a density (via the wf)
        //      returns the new F matrix

        //Get our new C's and D;s by diagonalization
        Eigen::SelfAdjointEigenSolver<matrix_type> eC(X.transpose()*F*X);
        matrix_type C = X*eC.eigenvectors().leftCols(nocc);
        D = C*C.transpose();

        //This really is only here for progress monitoring, in a scheme where
        //one iteration is a call to the
        std::cout<<"Iter: "<<iter
                 <<" dE: "<<dE
                 <<" [D,F]: "<<rmsComm
                 <<" E:"<<newE
                 <<std::endl;

        //Update in preperation for either returning or another cycle
        oldE=newE;
        newwfn=update_wfn(newwfn,&C,&D,&occs);

        //Check for convergence
        if(std::fabs(dE)<e_conv && std::fabs(rmsComm)<comm_conv)
            break;

        //Didn't converge, so next iteration
        ++iter;
    }while(iter<max_iters);

    //Returns our converged wavefunction and the corresponding energy
    return {newwfn,{oldE}};
}

}
