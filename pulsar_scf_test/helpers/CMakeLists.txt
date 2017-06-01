foreach(x TestDIIS TestShellPairItr TestShellTripleItr TestShellQuartetItr)
    pulsar_cxx_test(helpers ${x} ${PULSAR_SCF_PREFIX})
endforeach()

pulsar_py_test(helpers TestSchwarzScreen ${PULSAR_SCF_PREFIX})
