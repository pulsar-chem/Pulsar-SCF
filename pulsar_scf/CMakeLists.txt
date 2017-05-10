cmake_minimum_required(VERSION 3.1)
project(pulsar-scf CXX)
find_package(pulsar REQUIRED)
include(ExternalProject)

get_filename_component(LIBSCF_ROOT ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

add_library(pulsar_scf MODULE creator.cpp
                              D.cpp
                              DF.cpp
                              F.cpp
                              G.cpp
                              H.cpp
                              HelperFunctions.cpp
                              JK.cpp
                              S.cpp
                              SCF.cpp
                              SchwarzScreen.cpp
                              T.cpp
                              V.cpp
                              X.cpp
)
target_compile_options(pulsar_scf PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

target_include_directories(pulsar_scf PRIVATE pulsar ${LIBSCF_ROOT})
target_link_libraries(pulsar_scf PRIVATE pulsar)
install(TARGETS pulsar_scf
        LIBRARY DESTINATION ${PULSAR_MODULE_DIR}/pulsar_scf
)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        DESTINATION ${PULSAR_MODULE_DIR}
        FILES_MATCHING PATTERN "*.py"
)
