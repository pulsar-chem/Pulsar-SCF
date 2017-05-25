#include <pulsar_scf/helpers/ShellPairItr.hpp>
#include <pulsar/testing/CppTester.hpp>
#include "../TestCommon.hpp"
using namespace pulsar;
using namespace pulsar_scf;

std::vector<std::array<size_t,2>> corr_shells={
    {0,0},{1,0},{1,1},{2,0},{2,1},{2,2},{3,0},{3,1},{3,2},{3,3},{4,0},{4,1},
    {4,2},{4,3},{4,4}
};

std::map<std::array<size_t,2>,std::vector<std::array<size_t,2>>> corr_bfs={
  {{0,0},{{0,0}}},
  {{1,0},{{1,0}}},
  {{1,1},{{1,1}}},
  {{2,0},{{2,0},{3,0},{4,0}}},
  {{2,1},{{2,1},{3,1},{4,1}}},
  {{2,2},{{2,2},{2,3},{2,4},{3,2},{3,3},{3,4},{4,2},{4,3},{4,4}}},
  {{3,0},{{5,0}}},
  {{3,1},{{5,1}}},
  {{3,2},{{5,2},{5,3},{5,4}}},
  {{3,3},{{5,5}}},
  {{4,0},{{6,0}}},
  {{4,1},{{6,1}}},
  {{4,2},{{6,2},{6,3},{6,4}}},
  {{4,3},{{6,5}}},
  {{4,4},{{6,6}}}
};

TEST_SIMPLE(TestShellPairItr){
    CppTester tester("Testing the ShellPairItr classes");

    auto wf=make_wf();
    auto bs=wf.system->get_basis_set("PRIMARY");
    ShellPairItr my_itr(bs);
    auto corr_shell_itr=corr_shells.begin();
    while(my_itr)
    {
        const auto& idx=*my_itr;
        std::stringstream ss;
        ss<<"Shell pair: ( "<<(*corr_shell_itr)[0]<<" | "
          <<(*corr_shell_itr)[1]<<" )";
        tester.test_equal(ss.str(),idx,*corr_shell_itr);
        auto corr_bf_itr=corr_bfs.at(*corr_shell_itr).cbegin();
        for(const auto& bf_pair:my_itr)
        {
            std::stringstream ss1;
            ss1<<"Basis function pair: ( "<<(*corr_shell_itr)[0]<<": "
               <<(*corr_bf_itr)[0]<< " | "<<(*corr_shell_itr)[1]<<": "
               <<(*corr_bf_itr)[1]<<" )";
            tester.test_equal(ss1.str(),bf_pair,*corr_bf_itr);
            ++corr_bf_itr;
        }
        ++my_itr;++corr_shell_itr;
    }
    tester.print_results();
    return tester.nfailed();
}
