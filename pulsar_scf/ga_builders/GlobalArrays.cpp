#include "pulsar_scf/ga_builders/GlobalArrays.hpp"
#include <pulsar/util/Pybind11.hpp>
#include <ga.h>
#include <macdecls.h>

void* replace_malloc(size_t bytes, int align, char *name)
{
return malloc(bytes);
}

void replace_free(void *ptr)
{
free(ptr);
}

namespace pulsar_scf {
//Exists basically to move initialization/finalization to RAII model
struct GAEnv {
    GAEnv()
    {
        NGA_Initialize();
    }

    void MA_Init(int heap=-1,int stack=-1)
    {
        MA_init(C_DBL,stack,heap);
        GA_Register_stack_memory(replace_malloc, replace_free);
    }

    ~GAEnv()
    {
        //GA_Terminate();
    }
};

GAEnv Env_;

void GA_initialize(int heap,int stack){
    Env_.MA_Init(heap,stack);
}
void GA_finalize()
{
}

}//end namespace pulsar_scf

PYBIND11_PLUGIN(pulsar_scf)
{
    pybind11::module m("pulsar_scf", "PulsarChem SCF");
    m.def("GA_initialize", pulsar_scf::GA_initialize,
          pybind11::arg("heap") = -1,
          pybind11::arg("stack") = -1);
    m.def("GA_finalize",pulsar_scf::GA_finalize);
    return m.ptr();
}
