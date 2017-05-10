#pragma once
#include<pulsar/math/TensorImpl.hpp>
#include<vector>
#include<memory>

namespace pulsar_scf {
namespace detail_ {

template<typename matrix_type,typename int_ptr>
matrix_type fill_symmetric(const int_ptr& Ints,
                           const pulsar::BasisSet& bs)
{
    matrix_type result(bs.n_functions(),bs.n_functions());
    const size_t nshells1=bs.n_shell(),nshells2=bs.n_shell();
    for(size_t shell_i=0; shell_i<nshells1;++shell_i)
    {
        const size_t nbf_i=bs.shell(shell_i).n_functions();
        const size_t ni_off=bs.shell_start(shell_i);
        for(size_t shell_j=shell_i; shell_j<nshells2;++shell_j)
        {
            const size_t nbf_j=bs.shell(shell_j).n_functions();
            const size_t nj_off=bs.shell_start(shell_j);
            //Buffer is nbfi by nbfj
            const double* buffer=Ints->calculate(shell_i,shell_j);
            for(size_t mu=0;mu<nbf_i;++mu)
            {
                const size_t row_big=mu+ni_off;
                const size_t row_small=mu*nbf_j;
                for(size_t nu=0;nu<nbf_j;++nu)
                {
                    result(row_big,nu+nj_off)=buffer[row_small+nu];
                    result(nu+nj_off,row_big)=buffer[row_small+nu];
                }
            }
        }
    }
    return result;
}

//Untested
template<typename matrix_type,typename int_ptr>
matrix_type fill_asymmetric(const int_ptr& Ints,
                            const pulsar::BasisSet& bs1,
                            const pulsar::BasisSet& bs2)
{
    matrix_type result(bs1.n_functions(),bs2.n_functions());
    const size_t nshells1=bs1.n_shell(),nshells2=bs2.n_shell();
    for(size_t shell_i=0; shell_i<nshells1;++shell_i)
    {
        const size_t nbf_i=bs1.shell(shell_i).n_functions();
        const size_t ni_off=bs1.shell_start(shell_i);
        for(size_t shell_j=0; shell_j<nshells2;++shell_j)
        {
            const size_t nbf_j=bs2.shell(shell_j).n_functions();
            const size_t nj_off=bs2.shell_start(shell_j);
            //Buffer is nbfi by nbfj
            const double* buffer=Ints->calculate(shell_i,shell_j);
            for(size_t mu=0;mu<nbf_i;++mu)
            {
                const size_t row_big=mu+ni_off;
                const size_t row_small=mu*nbf_j;
                for(size_t nu=0;nu<nbf_j;++nu)
                    result(row_big,nu+nj_off)=buffer[row_small+nu];
            }
        }
    }
    return result;

}


}}//End namespaces
