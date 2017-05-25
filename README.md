# Pulsar-SCF

This is the beginning of a prototyping effort for NWChemEX's SCF code.  As it
currently stands this code is heavily borrowed from the example in LibInt's test
directory.  The goal is over the next several months to significantly rewrite it
so that is feature rich, more optimized, and parallelized.  I also intend to try
several tensor libraries out so look for tensor library specific implementations
to come.

## Building

Pulsar-SCF uses a CMake superbuild.  For the most part, this means the following
will be sufficient to build Pulsar-SCF:

~~~.sh
cmake -H. -Bbuild -DCMAKE_VARIABLE1=value -DCMAKE_VARIABLE2=...
cd build && make
make install
~~~

--------------------------------------------------------------------------------
| CMake Variable | Description                                                 |
| :------------: | :-----------------------------------------------------------|
| CMAKE_CXX_COMPILER | The C++ compiler that will be used                      |
| CMAKE_CXX_FLAGS | Flags that will be passed to C++ compiler                  |
| CMAKE_C_COMPILER | The C compiler that will be used                          |
| CMAKE_C_FLAGS | Flags that will be passed to the C compiler                  |
| MPI_C_COMPILER | The MPI wrapper compiler for building C executables         |
| MPI_CXX_COMPILER | The MPI wrapper compiler for building C++ executables     |
| CMAKE_BUILD_TYPE | Debug, Release, or RelWithDebInfo                         |
| CMAKE_PREFIX_PATH | A list of places CMake will look for dependencies        |
| CMAKE_INSTALL_PREFIX | The install directory                                 |
--------------------------------------------------------------------------------

By default, the only dependency of Pulsar-SCF is Pulsar-Core.  If you enable
additional features more dependencies may be introduced.  Consult the next
section for more details.

### Optional Features

The following is a table of all optional features of Pulsar-SCF.  Following the
table is a more thorough description of each option.

-------------------------------------------------------------------
| Option Name | Default | Description                             |
|     :---:   |  :---:  | :---------------------------------------|
| ENABLE_GA   | False   | Should global array support be enabled? |
-------------------------------------------------------------------

- `ENABLE_GA` This controls whether a version of the SCF is built that uses the
  Global Arrays.  If enabled, the build will attempt to locate an existing copy
  of Global Arrays and if it does not find one will build one.  Note, Global
  arrays must have been built with `-fPIC` and `--enable-cxx` to work with
  Pulsar.

## Documentation

TODO: Add when automated
