foreach(x TestD TestDFCoef TestDFInts TestF TestG TestH TestJK TestMetric TestS TestT TestV TestX)
   pulsar_py_test(builders ${x} ${PULSAR_SCF})
endforeach()
