#include "pulsar_scf/ga_builders/GlobalArrays.hpp"
#include "pulsar_scf/ga_builders/GATensor.hpp"
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
       //TODO: Set MPI comm via this function
        NGA_Initialize();
    }

    void MA_Init(int heap=-1,int stack=-1)
    {
        MA_init(C_DBL,stack,heap);
        GA_Register_stack_memory(replace_malloc, replace_free);
    }

    ~GAEnv()
    {
                GA_Terminate();
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

//Function for exporting our GA adapted tensors
template<size_t N>
void export_GA_x_impl(pybind11::module& m,const char* Name)
{
    using namespace pulsar_scf;
    using GA_t=GATensorImpl<N>;
    pybind11::class_<GA_t,pulsar::TensorImpl<N,double>,std::shared_ptr<GA_t>>(m,Name)
    .def(pybind11::init<const GA_t&>())
    .def(pybind11::self == pybind11::self)
    .def(pybind11::self != pybind11::self)
    .def("my_hash",&GA_t::my_hash)
    .def("sizes",&GA_t::sizes)
    .def("get_value",&GA_t::get_value)
    .def("set_value",&GA_t::set_value)
    ;
}


PYBIND11_PLUGIN(pulsar_scf)
{
    pybind11::module m("pulsar_scf", "PulsarChem SCF");
    m.def("GA_initialize", pulsar_scf::GA_initialize,
          pybind11::arg("heap") = -1,
          pybind11::arg("stack") = -1);
    m.def("GA_finalize",pulsar_scf::GA_finalize);

    export_GA_x_impl<2>(m,"GAMatrixImpl");
    export_GA_x_impl<3>(m,"GARank3Impl");

    return m.ptr();
}
