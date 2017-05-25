include(ExternalProject)

if(ENABLE_GA)
    find_package(GlobalArrays)
    find_package(MPI)
    if(NOT GLOBALARRAYS_FOUND)
        message(STATUS "GA NOT FOUND")
        set(GA_PREFIX ${CMAKE_BINARY_DIR}/ga)
        set(GA_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        ExternalProject_Add(global_arrays
            PREFIX ${GA_PREFIX}
            GIT_REPOSITORY https://github.com/GlobalArrays/ga.git
            GIT_TAG v5.6
            UPDATE_COMMAND ./autogen.sh
            CONFIGURE_COMMAND ${GA_PREFIX}/src/global_arrays/configure
                              CXX=${MPI_CXX_COMPILER}
                              CC=${MPI_C_COMPILER}
                              CXXFLAGS=${GA_FLAGS}
                              CFLAGS=${GA_FLAGS}
                              --with-pic
                              --prefix=${STAGE_INSTALL_PREFIX}
                              --enable-cxx
            BUILD_COMMAND make
            INSTALL_DIR ${STAGE_INSTALL_PREFIX}
            INSTALL_COMMAND make install
        )
    else()#Already built one found
        add_library(global_arrays INTERFACE)
    endif()
else()#Not building
   add_library(global_arrays INTERFACE)
endif()
