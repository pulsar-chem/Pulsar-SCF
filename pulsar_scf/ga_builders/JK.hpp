#pragma once
#include <pulsar/modulebase/MatrixBuilder.hpp>
#include "pulsar_scf/ga_builders/GlobalArrays.hpp"
#include "pulsar_scf/PulsarSCF.hpp"

namespace pulsar_scf {

#ifdef ENABLE_GA

MATRIX_BUILDER(GAJK)

#else

NO_GA_MATRIX_BUILDER(GAJK)

#endif


}//End namespace pulsar_scf
