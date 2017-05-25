
if(${ENABLE_GA})
    set(PULSAR_GA_BUILDERS_FILES ${PULSAR_GA_BUILDERS_FILES}
                                 GATensor.cpp
                                 GlobalArrays.cpp
                                 dummy.F
                                 JK.cpp
                                 PARENT_SCOPE
    )
    set_source_files_properties(pulsar_scf/ga_builders/JK.cpp PROPERTIES COMPILE_FLAGS -DENABLE_GA)
endif()
