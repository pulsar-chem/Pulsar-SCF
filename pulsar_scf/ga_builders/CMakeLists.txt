
if(${ENABLE_GA})
    set(PULSAR_GA_BUILDERS_FILES ${PULSAR_GA_BUILDERS_FILES}
                                 GATensor.cpp
                                 GlobalArrays.cpp
                                 dummy.F
                                 H.cpp
                                 G.cpp
                                 F.cpp
                                 JK.cpp
                                 T.cpp
                                 V.cpp
                                 S.cpp
                                 X.cpp
                                 PARENT_SCOPE
    )
    set_source_files_properties(pulsar_scf/ga_builders/JK.cpp PROPERTIES COMPILE_FLAGS -DENABLE_GA)
endif()
