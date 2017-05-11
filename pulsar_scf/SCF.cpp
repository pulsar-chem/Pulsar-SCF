#include "pulsar_scf/SCF.hpp"
#include "pulsar_scf/helpers/DIIS.hpp"
#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/modulebase/MatrixBuilder.hpp>

using namespace pulsar;
using matrix_type=EigenMatrixImpl::matrix_type;

double nuclear_nuclear_repulsion(const Wavefunction& wfn)
{
    double ENuc=0.0;
    for(auto ai: *wfn.system)
    {
        for(auto aj: *wfn.system)
        {
            if(aj==ai)break;
            const double dx=ai[0]-aj[0],dy=ai[1]-aj[1],dz=ai[2]-aj[2];
            ENuc+=ai.Z*aj.Z/std::sqrt(dx*dx+dy*dy+dz*dz);
        }
    }
    return ENuc;
}

namespace pulsar_scf {

DerivReturnType SCF::deriv_(size_t deriv,
                            const Wavefunction & wfn)
{
    auto bs=wfn.system->get_basis_set("PRIMARY");
    const double e_conv=1e-6;
    const double comm_conv=1e-6;
    const size_t max_iters=100;
    size_t iter=0;

    const size_t nocc=5;

    const double ENuc=nuclear_nuclear_repulsion(wfn);
    double oldE=ENuc;
    std::cout<<"Nuclear-Nuclear repulsion is: "<<ENuc<<std::endl;
    std::unique_ptr<DIIS<matrix_type>> diis;
    const auto H=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("H_KEY")
                                   ->calculate("",deriv,wfn,bs,bs)[0]);
    const auto S=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("S_KEY")
                                   ->calculate("",deriv,wfn,bs,bs)[0]);
    const auto Xs=create_child_from_option<MatrixBuilder>("X_KEY")->calculate("",deriv,wfn,bs,bs);
    matrix_type X=*convert_to_eigen(*Xs[0]);
    matrix_type D=*convert_to_eigen(*wfn.opdm->get(Irrep::A,Spin::alpha));
    Wavefunction newwfn(wfn);
    do{
        auto F=*convert_to_eigen(*create_child_from_option<MatrixBuilder>("F_KEY")
                     ->calculate("",deriv,newwfn,bs,bs)[0]);
        const double newE=D.cwiseProduct(H+F).sum()+ENuc;
        const double dE=newE-oldE;
        const auto comm=(F*D*S-S*D*F);
        const double rmsComm=comm.norm();
        if(!diis)//First iteration
            diis=std::unique_ptr<DIIS<matrix_type>>(new DIIS<matrix_type>(F,comm));
        else F=diis->update(F,comm);
        Eigen::SelfAdjointEigenSolver<matrix_type> eC(X.transpose()*F*X);
        matrix_type C = X*eC.eigenvectors().leftCols(nocc);
        D = C*C.transpose();
        std::cout<<"Iter: "<<iter
                 <<" dE: "<<dE
                 <<" [D,F]: "<<rmsComm
                 <<" E:"<<newE
                 <<std::endl;

        if(std::fabs(dE)<e_conv && std::fabs(rmsComm)<comm_conv)break;
        oldE=newE;
        newwfn=update_wfn(newwfn,&C,&D);
        ++iter;
    }while(iter<max_iters);
    return {newwfn,{oldE}};
}

}
