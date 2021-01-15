#include <catch.hpp>

#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <mockturtle/algorithms/exact_syn/sat.hpp>

TEST_CASE( "Synthesis with SAT using partial DAG", "[aqfp_network_sat]" )
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  Ntk pdag1 = { { 3, 3, 0, 0, 0, 0, 0 }, { { 1, 2, 3 }, { 4, 5, 6 }, {}, {}, {}, {}, {} }, { 2, 3, 4, 5, 6 }, 0, 0u, true, 2u, -1.0, { { 3, 2 } }, {}, {} };

  std::vector<uint64_t> tts = { 0x01ul, 0xe8ul };
  for ( auto tt : tts )
  {
    mockturtle::sat_encoder se( pdag1, 3u, tt, 0u );
    CHECK( se.solve() );
  }

  mockturtle::sat_encoder se( pdag1, 3u, 0x69ul, 0u );
  CHECK( !se.solve() );
}

TEST_CASE( "Synthesis with SAT using DAG", "[aqfp_network_sat]" )
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  Ntk dag1 = { { 3, 3, 0, 0, 0}, { { 1, 2, 3 }, { 2, 3, 4 }, {}, {}, {} }, { 2, 3, 4}, 2, 0u, false, 2u, -1.0, { { 3, 2 } }, {}, {} };

  std::vector<uint64_t> tts = { 0x69ul, 0xe8ul };
  for ( auto tt : tts )
  {
    mockturtle::sat_encoder se( dag1, 3u, tt, 0u );
    CHECK( !se.solve() );
  }

  mockturtle::sat_encoder se( dag1, 2u, 0x1ul, 0u );
  CHECK( se.solve() );
}