#include <catch.hpp>

#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <mockturtle/algorithms/exact_syn/simulate_dag.hpp>

TEST_CASE( "Simulate AQFP logical networks", "[aqfp_network_sim]" )
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  mockturtle::dag_simulator<uint64_t> sim( { 0x00ul, 0xaaul, 0xccul, 0xf0ul } );

  Ntk dag1 = { { 3, 3, 0, 0, 0 }, { { 1, 2, 3 }, { 2, 3, 4 }, {}, {}, {} }, { 2, 3, 4 }, 0, 0u, 0, 2u, -1, {  {3, 2 } }, {}, {} };
  auto t1 = sim.all_functions_from_dag( dag1 );
  std::unordered_set<uint64_t> s1 = { 142ul, 204ul, 178ul, 232ul, 212ul, 170ul };
  CHECK( s1 == t1 );

  Ntk dag2 = { { 3, 3, 0, 0, 0 }, { { 1, 2, 3 }, { 2, 3, 4 }, {}, {}, {} }, { 2, 3, 4 }, 3, 0u, 0, 2u, -1, {  { 3, 2 } }, {}, {} };
  auto t2 = sim.all_functions_from_dag( dag2 );
  std::unordered_set<uint64_t> s2 = { 34, 136, 238, 68, 0, 170 };

  CHECK( s2 == t2 );
}
