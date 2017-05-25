cmake_minimum_required(VERSION 3.1)
project(pulsar-scf CXX)
find_package(pulsar REQUIRED)
include(ExternalProject)

get_filename_component(LIBSCF_ROOT ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
foreach(sub_dir builders ga_builders helpers)

    # Puts all the files into a list called PULSAR_${SUB_DIR}_FILES
    add_subdirectory(${sub_dir})

    # Make the list a full path and add to the master list
    string(TOUPPER ${sub_dir} SUB_DIR)
    foreach(F ${PULSAR_${SUB_DIR}_FILES})
        list(APPEND PULSAR_SCF_FILE_LIST "${sub_dir}/${F}")
    endforeach()
endforeach()

add_library(pulsar_scf MODULE creator.cpp
                              HelperFunctions.cpp
                              SCF.cpp
                              ${PULSAR_SCF_FILE_LIST}
)
target_compile_options(pulsar_scf PRIVATE ${PULSAR_CXX_STRICT_FLAGS})

target_include_directories(pulsar_scf PRIVATE pulsar ${LIBSCF_ROOT})
target_link_libraries(pulsar_scf PRIVATE pulsar)

if(ENABLE_GA)
    enable_language(Fortran)
    find_package(GlobalArrays REQUIRED)
    target_compile_definitions(pulsar_scf PRIVATE ENABLE_GA)
    target_include_directories(pulsar_scf PRIVATE ${GLOBALARRAYS_INCLUDE_DIRS})
    target_link_libraries(pulsar_scf PRIVATE ${GLOBALARRAYS_LIBRARIES})
endif()

install(TARGETS pulsar_scf
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/pulsar_scf
)
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
        DESTINATION ${CMAKE_INSTALL_PREFIX}
        FILES_MATCHING PATTERN "*.py"
)
