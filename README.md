# Pulsar-SCF

This is the beginning of a prototyping effort for NWChemEX's SCF code.  As it currently
stands this code is heavily borrowed from the example in LibInt's test directory.  The
goal is over the next several months to significantly rewrite it so that is feature
rich, more optimized, and parallelized.  I also intend to try several tensor libraries
out so look for tensor library specific implementations to come.

## Building

Whatever CMake command you ran to build Pulsar-Core should be sufficient to configure Pulsar-SCF
as it has no dependencies aside from Pulsar-Core (and Pulsar-Core's dependencies).
