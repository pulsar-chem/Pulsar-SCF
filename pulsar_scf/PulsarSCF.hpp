#pragma once

///A macro for declaring a new MatrixBuilder module
#define MATRIX_BUILDER(name)\
    class name: public pulsar::MatrixBuilder {\
    public:\
        name(ID_t id):MatrixBuilder(id){}\
        ReturnType calculate_(const std::string &,\
                              unsigned int deriv,\
                              const pulsar::Wavefunction &,\
                              const pulsar::BasisSet &,\
                              const pulsar::BasisSet &);\
        HashType my_hash_(const std::string &,\
                           unsigned int,\
                          const pulsar::Wavefunction &,\
                          const pulsar::BasisSet &,\
                          const pulsar::BasisSet &);\
    };

///A macro for declaring a new Rank3Builder
#define RANK3_BUILDER(name)\
    class name: public pulsar::Rank3Builder {\
    public:\
        name(ID_t id):Rank3Builder(id){}\
        ReturnType calculate_(const std::string & key,\
                              unsigned int deriv,\
                              const pulsar::Wavefunction & wfn,\
                              const pulsar::BasisSet & bs1,\
                              const pulsar::BasisSet & bs2,\
                              const pulsar::BasisSet & bs3);\
    \
        HashType my_hash_(const std::string & key,\
                          unsigned int deriv,\
                          const pulsar::Wavefunction &wfn,\
                          const pulsar::BasisSet &bs1,\
                          const pulsar::BasisSet &bs2,\
                          const pulsar::BasisSet &bs3);\
    };
