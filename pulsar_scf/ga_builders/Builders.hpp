#pragma once
#include "pulsar_scf/ga_builders/GlobalArrays.hpp"

namespace pulsar_scf {

///Builds J and K via a direct algorithm
GA_MATRIX_BUILDER(GAJK)

///Builds G
GA_MATRIX_BUILDER(GAG)

///Builds F
GA_MATRIX_BUILDER(GAF)

///Builds the kinetic energy of the electrons integrals
GA_MATRIX_BUILDER(GAT)

///Builds the nucleus electron potential integrals
GA_MATRIX_BUILDER(GAV)

///Builds the overlap matrix
GA_MATRIX_BUILDER(GAS)

///Builds the orthogonalizer matrix
GA_MATRIX_BUILDER(GAX)

///Builds the core Hamiltonian
GA_MATRIX_BUILDER(GAH)

}//End namespace pulsar_scf
