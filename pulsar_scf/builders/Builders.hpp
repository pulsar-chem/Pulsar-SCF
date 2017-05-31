#pragma once
#include "pulsar_scf/PulsarSCF.hpp"

namespace pulsar_scf {

///Builds the overlap matrix
MATRIX_BUILDER(Overlap)

///Builds the orthogonalizer for the Fock Matrix
MATRIX_BUILDER(Orthogonalizer)

///Builds the Fock matrix
MATRIX_BUILDER(F)

///Builds the 2e part of the Fock Matrix
MATRIX_BUILDER(G)

///Builds the J and K matrices via a direct algorithm
MATRIX_BUILDER(JK)

///Builds the J and K matrices via a density-fit core algorithm
MATRIX_BUILDER(DFJK)

///Builds the 1e part of the Fock Matrix
MATRIX_BUILDER(HCore)

///Builds the kinetic energy of the electrons integrals
MATRIX_BUILDER(TElectronic)

///Builds the nucleus electron attraction integrals
MATRIX_BUILDER(NuclearElectronic)

///Builds \f$[J]_{PQ}\f$
MATRIX_BUILDER(Metric)

///Builds integrals (mu nu | Q)
RANK3_BUILDER(DFInts)

///Builds d_Qmn
RANK3_BUILDER(DFCoef)

}
