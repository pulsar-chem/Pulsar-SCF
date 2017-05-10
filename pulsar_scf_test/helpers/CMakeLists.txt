pulsar_cxx_test(helpers TestDIIS ${PULSAR_SCF})
target_link_libraries(TestDIIS PRIVATE /home/ryan/pulsar/lib/libint2.a)
pulsar_py_test(helpers TestSchwarzScreen ${PULSAR_SCF})
pulsar_cxx_test(helpers TestShellPairItr ${PULSAR_SCF})
