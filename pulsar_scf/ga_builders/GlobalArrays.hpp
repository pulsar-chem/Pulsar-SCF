#pragma once
#include<pulsar/exception/PulsarException.hpp>
#include"pulsar_scf/PulsarSCF.hpp"

namespace pulsar_scf {

///Initializes GA with the given heap and stack sizes (negative means GA picks)
void GA_initialize(int heap=-1, int stack=-1);
///Finalizes GA
void GA_finalize();

}

#ifdef ENABLE_GA
#define GA_MATRIX_BUILDER(name) MATRIX_BUILDER(name)
#else
#define GA_MATRIX_BUILDER(name)\
    class name: public pulsar::MatrixBuilder {\
    public:\
        name(ID_t id):MatrixBuilder(id){}\
        ReturnType calculate_(const std::string &,\
                              unsigned int deriv,\
                              const pulsar::Wavefunction &,\
                              const pulsar::BasisSet &,\
                              const pulsar::BasisSet &){\
         throw pulsar::PulsarException("GA was not enabled.");\
        }\
        HashType my_hash_(const std::string &,\
                           unsigned int,\
                          const pulsar::Wavefunction &,\
                          const pulsar::BasisSet &,\
                          const pulsar::BasisSet &){\
          throw pulsar::PulsarException("GA was not enabled.");\
        }\
    };
#endif
