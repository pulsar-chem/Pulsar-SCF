#include <pulsar_scf/helpers/ShellTripleItr.hpp>
#include <pulsar/testing/CppTester.hpp>
#include "../TestCommon.hpp"
using namespace pulsar;
using namespace pulsar_scf;

std::vector<std::array<size_t,3>> corr_shells={
    {0,0,0},{0,1,0},{0,1,1},{0,2,0},{0,2,1},{0,2,2},{0,3,0},{0,3,1},{0,3,2},
    {0,3,3},{0,4,0},{0,4,1},{0,4,2},{0,4,3},{0,4,4},{1,0,0},{1,1,0},{1,1,1},
    {1,2,0},{1,2,1},{1,2,2},{1,3,0},{1,3,1},{1,3,2},{1,3,3},{1,4,0},{1,4,1},
    {1,4,2},{1,4,3},{1,4,4},{2,0,0},{2,1,0},{2,1,1},{2,2,0},{2,2,1},{2,2,2},
    {2,3,0},{2,3,1},{2,3,2},{2,3,3},{2,4,0},{2,4,1},{2,4,2},{2,4,3},{2,4,4},
    {3,0,0},{3,1,0},{3,1,1},{3,2,0},{3,2,1},{3,2,2},{3,3,0},{3,3,1},{3,3,2},
    {3,3,3},{3,4,0},{3,4,1},{3,4,2},{3,4,3},{3,4,4},{4,0,0},{4,1,0},{4,1,1},
    {4,2,0},{4,2,1},{4,2,2},{4,3,0},{4,3,1},{4,3,2},{4,3,3},{4,4,0},{4,4,1},
    {4,4,2},{4,4,3},{4,4,4},
};

std::map<std::array<size_t,3>,std::vector<std::array<size_t,3>>> corr_bfs={
    {{0,0,0},{{0,0,0},}},
    {{0,1,0},{{0,1,0},}},
    {{0,1,1},{{0,1,1},}},
    {{0,2,0},{{0,2,0},{0,3,0},{0,4,0},}},
    {{0,2,1},{{0,2,1},{0,3,1},{0,4,1},}},
    {{0,2,2},{{0,2,2},{0,2,3},{0,2,4},{0,3,2},{0,3,3},{0,3,4},{0,4,2},{0,4,3},{0,4,4},}},
    {{0,3,0},{{0,5,0},}},
    {{0,3,1},{{0,5,1},}},
    {{0,3,2},{{0,5,2},{0,5,3},{0,5,4},}},
    {{0,3,3},{{0,5,5},}},
    {{0,4,0},{{0,6,0},}},
    {{0,4,1},{{0,6,1},}},
    {{0,4,2},{{0,6,2},{0,6,3},{0,6,4},}},
    {{0,4,3},{{0,6,5},}},
    {{0,4,4},{{0,6,6},}},
    {{1,0,0},{{1,0,0},}},
    {{1,1,0},{{1,1,0},}},
    {{1,1,1},{{1,1,1},}},
    {{1,2,0},{{1,2,0},{1,3,0},{1,4,0},}},
    {{1,2,1},{{1,2,1},{1,3,1},{1,4,1},}},
    {{1,2,2},{{1,2,2},{1,2,3},{1,2,4},{1,3,2},{1,3,3},{1,3,4},{1,4,2},{1,4,3},{1,4,4},}},
    {{1,3,0},{{1,5,0},}},
    {{1,3,1},{{1,5,1},}},
    {{1,3,2},{{1,5,2},{1,5,3},{1,5,4},}},
    {{1,3,3},{{1,5,5},}},
    {{1,4,0},{{1,6,0},}},
    {{1,4,1},{{1,6,1},}},
    {{1,4,2},{{1,6,2},{1,6,3},{1,6,4},}},
    {{1,4,3},{{1,6,5},}},
    {{1,4,4},{{1,6,6},}},
    {{2,0,0},{{2,0,0},{3,0,0},{4,0,0},}},
    {{2,1,0},{{2,1,0},{3,1,0},{4,1,0},}},
    {{2,1,1},{{2,1,1},{3,1,1},{4,1,1},}},
    {{2,2,0},{{2,2,0},{2,3,0},{2,4,0},{3,2,0},{3,3,0},{3,4,0},{4,2,0},{4,3,0},{4,4,0},}},
    {{2,2,1},{{2,2,1},{2,3,1},{2,4,1},{3,2,1},{3,3,1},{3,4,1},{4,2,1},{4,3,1},{4,4,1},}},
    {{2,2,2},{{2,2,2},{2,2,3},{2,2,4},{2,3,2},{2,3,3},{2,3,4},{2,4,2},{2,4,3},{2,4,4},{3,2,2},{3,2,3},{3,2,4},{3,3,2},{3,3,3},{3,3,4},{3,4,2},{3,4,3},{3,4,4},{4,2,2},{4,2,3},{4,2,4},{4,3,2},{4,3,3},{4,3,4},{4,4,2},{4,4,3},{4,4,4},}},
    {{2,3,0},{{2,5,0},{3,5,0},{4,5,0},}},
    {{2,3,1},{{2,5,1},{3,5,1},{4,5,1},}},
    {{2,3,2},{{2,5,2},{2,5,3},{2,5,4},{3,5,2},{3,5,3},{3,5,4},{4,5,2},{4,5,3},{4,5,4},}},
    {{2,3,3},{{2,5,5},{3,5,5},{4,5,5},}},
    {{2,4,0},{{2,6,0},{3,6,0},{4,6,0},}},
    {{2,4,1},{{2,6,1},{3,6,1},{4,6,1},}},
    {{2,4,2},{{2,6,2},{2,6,3},{2,6,4},{3,6,2},{3,6,3},{3,6,4},{4,6,2},{4,6,3},{4,6,4},}},
    {{2,4,3},{{2,6,5},{3,6,5},{4,6,5},}},
    {{2,4,4},{{2,6,6},{3,6,6},{4,6,6},}},
    {{3,0,0},{{5,0,0},}},
    {{3,1,0},{{5,1,0},}},
    {{3,1,1},{{5,1,1},}},
    {{3,2,0},{{5,2,0},{5,3,0},{5,4,0},}},
    {{3,2,1},{{5,2,1},{5,3,1},{5,4,1},}},
    {{3,2,2},{{5,2,2},{5,2,3},{5,2,4},{5,3,2},{5,3,3},{5,3,4},{5,4,2},{5,4,3},{5,4,4},}},
    {{3,3,0},{{5,5,0},}},
    {{3,3,1},{{5,5,1},}},
    {{3,3,2},{{5,5,2},{5,5,3},{5,5,4},}},
    {{3,3,3},{{5,5,5},}},
    {{3,4,0},{{5,6,0},}},
    {{3,4,1},{{5,6,1},}},
    {{3,4,2},{{5,6,2},{5,6,3},{5,6,4},}},
    {{3,4,3},{{5,6,5},}},
    {{3,4,4},{{5,6,6},}},
    {{4,0,0},{{6,0,0},}},
    {{4,1,0},{{6,1,0},}},
    {{4,1,1},{{6,1,1},}},
    {{4,2,0},{{6,2,0},{6,3,0},{6,4,0},}},
    {{4,2,1},{{6,2,1},{6,3,1},{6,4,1},}},
    {{4,2,2},{{6,2,2},{6,2,3},{6,2,4},{6,3,2},{6,3,3},{6,3,4},{6,4,2},{6,4,3},{6,4,4},}},
    {{4,3,0},{{6,5,0},}},
    {{4,3,1},{{6,5,1},}},
    {{4,3,2},{{6,5,2},{6,5,3},{6,5,4},}},
    {{4,3,3},{{6,5,5},}},
    {{4,4,0},{{6,6,0},}},
    {{4,4,1},{{6,6,1},}},
    {{4,4,2},{{6,6,2},{6,6,3},{6,6,4},}},
    {{4,4,3},{{6,6,5},}},
    {{4,4,4},{{6,6,6},}},
};

TEST_SIMPLE(TestShellTripleItr){
    CppTester tester("Testing the ShellTripleItr classes");

    auto wf=make_wf();
    auto bs=wf.system->get_basis_set("PRIMARY");
    ShellTripleItr my_itr(bs,bs);
    auto corr_shell_itr=corr_shells.begin();
    while(my_itr)
    {
        const auto& idx=*my_itr;
        std::stringstream ss;
        ss<<"Shell triple: ( "<<(*corr_shell_itr)[0]<<" | "
          <<(*corr_shell_itr)[1]<<" "<<(*corr_shell_itr)[2]<<")";
        tester.test_equal(ss.str(),idx,*corr_shell_itr);
        auto corr_bf_itr=corr_bfs.at(*corr_shell_itr).cbegin();
        for(const auto& bf_pair:my_itr)
        {
            std::stringstream ss1;
            ss1<<"Basis function triple: ( "
               <<(*corr_shell_itr)[0]<<": "<<(*corr_bf_itr)[0]
               << " | "
               <<(*corr_shell_itr)[1]<<": "<<(*corr_bf_itr)[1]<<" "
               <<(*corr_shell_itr)[2]<<": "<<(*corr_bf_itr)[2]<<" )";
            tester.test_equal(ss1.str(),bf_pair,*corr_bf_itr);
            ++corr_bf_itr;
        }
        ++my_itr;++corr_shell_itr;
    }
    tester.print_results();
    return tester.nfailed();
}
