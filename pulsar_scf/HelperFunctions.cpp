#include<pulsar_scf/HelperFunctions.hpp>

using namespace pulsar;
using vector_type=EigenVectorImpl::vector_type;
using matrix_type=EigenMatrixImpl::matrix_type;


namespace pulsar_scf {

std::shared_ptr<const IrrepSpinVectorD> guess_occ(const Wavefunction& wf)
{
    if(wf.occupations)return wf.occupations;
    auto occs=std::make_shared<IrrepSpinVectorD>();
    const size_t nelectrons=wf.system->get_sum_n_electrons();
    const bool odd_electrons=nelectrons%2;
    if(odd_electrons && wf.system->multiplicity==1.0)
    {
        vector_type occa=vector_type::Zero(nelectrons);
        for(size_t i=0;i<nelectrons/2;++i)occa(i)=1.0;
        occs->set(Irrep::A,Spin::alpha,std::make_shared<EigenVectorImpl>(std::move(occa)));
    }
    return occs;
}

Wavefunction update_wfn(const Wavefunction& wf,
                        const matrix_type* Ca,
                        const matrix_type* Da,
                        const matrix_type* Cb,
                        const matrix_type* Db)
{
    Wavefunction newwfn(wf);
    IrrepSpinMatrixD opdm,cmat;
    if(Da)
        opdm.set(Irrep::A,Spin::alpha,std::make_shared<EigenMatrixImpl>(*Da));
    if(Db)
        opdm.set(Irrep::A,Spin::beta,std::make_shared<EigenMatrixImpl>(*Db));
    if(Ca)
        cmat.set(Irrep::A,Spin::alpha,std::make_shared<EigenMatrixImpl>(*Ca));
    if(Cb)
        cmat.set(Irrep::A,Spin::alpha,std::make_shared<EigenMatrixImpl>(*Cb));
    newwfn.opdm=std::make_shared<const IrrepSpinMatrixD>(std::move(opdm));
    newwfn.cmat=std::make_shared<const IrrepSpinMatrixD>(std::move(cmat));
    return newwfn;
}

}//End namespace
