#include "pulsar_scf/HelperFunctions.hpp"
#include <pulsar/exception/PulsarException.hpp>
#include "pulsar_scf/MatrixFillFxns.hpp"

using namespace pulsar;
using namespace std;

using vector_type=EigenVectorImpl::vector_type;
using matrix_type=EigenMatrixImpl::matrix_type;


namespace pulsar_scf {

size_t nshells(const BasisSet& bs)
{
    size_t nshells=0;
    for(const auto& shelli : bs)
        nshells+=shelli.n_general_contractions();
    return nshells;
}


shared_ptr<const IrrepSpinVectorD> guess_occ(const Wavefunction& wf)
{
    if(wf.occupations)return wf.occupations;
    auto occs=make_shared<IrrepSpinVectorD>();
    const size_t nelectrons=wf.system->get_sum_n_electrons();
    const bool odd_electrons=nelectrons%2;
    if(!odd_electrons && wf.system->multiplicity==1.0)
    {
        vector_type occa=vector_type::Zero(nelectrons);
        for(size_t i=0;i<nelectrons/2;++i)occa(i)=1.0;
        occs->set(Irrep::A,Spin::alpha,make_shared<EigenVectorImpl>(move(occa)));
    }
    else
    {
        throw PulsarException("I don't know how to guess the occupation for your system");
    }
    return occs;
}

Wavefunction update_wfn(const Wavefunction& wf,
                        const matrix_type* Ca,
                        const matrix_type* Da,
                        const vector_type* Epsa,
                        const matrix_type* Cb,
                        const matrix_type* Db,
                        const vector_type* Epsb)
{
    Wavefunction newwfn(wf);
    IrrepSpinMatrixD opdm,cmat;
    IrrepSpinVectorD epsilons;
    if(Da)
        opdm.set(Irrep::A,Spin::alpha,make_shared<EigenMatrixImpl>(*Da));
    if(Db)
        opdm.set(Irrep::A,Spin::beta,make_shared<EigenMatrixImpl>(*Db));
    if(Ca)
        cmat.set(Irrep::A,Spin::alpha,make_shared<EigenMatrixImpl>(*Ca));
    if(Cb)
        cmat.set(Irrep::A,Spin::beta,make_shared<EigenMatrixImpl>(*Cb));
    if(Epsa)
        epsilons.set(Irrep::A,Spin::alpha,make_shared<EigenVectorImpl>(*Epsa));
    if(Epsb)
        epsilons.set(Irrep::A,Spin::beta,make_shared<EigenVectorImpl>(*Epsb));
    newwfn.opdm=make_shared<const IrrepSpinMatrixD>(move(opdm));
    newwfn.cmat=make_shared<const IrrepSpinMatrixD>(move(cmat));
    newwfn.epsilon=make_shared<const IrrepSpinVectorD>(move(epsilons));
    return newwfn;
}

MatrixBuilder::ReturnType matrix_builder_kernel(unsigned int deriv,
                                                const Wavefunction& wfn,
                                                const BasisSet& bs1,
                                                const BasisSet& bs2,
                                                TwoCenterIntegral& Ints)
{
    Ints.initialize(deriv,wfn,bs1,bs2);
    const bool is_symmetric= bs1==bs2;
    auto result(move(is_symmetric ?
                detail_::fill_symmetric<matrix_type>(&Ints,bs1) :
                detail_::fill_asymmetric<matrix_type>(&Ints,bs1,bs2)
    ));
    return {make_shared<EigenMatrixImpl>(move(result))};
}

}//End namespace
