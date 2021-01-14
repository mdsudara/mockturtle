#include <catch.hpp>

#include <set>

#include <mockturtle/algorithms/exact_syn/dag_cost.hpp>

using namespace mockturtle;

mockturtle::aqfp_logical_network_t<int> net1 = {
    { 3, 3, 3, 3, 0, 0, 0, 0, 0 },
    { { 1, 4, 5 }, { 2, 4, 5 }, { 3, 5, 6 }, { 5, 7, 8 }, {}, {}, {}, {}, {} },
    { 4, 5, 6, 7, 8 },
    5u,

    0u,
    false,
    4u,
    -1.0,
    { { 3u, 4u } },
    {},
    {} };

mockturtle::aqfp_logical_network_t<int> net2 = {
    { 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0 },
    { { 1, 2, 6 }, { 4, 5, 6 }, { 3, 6, 7 }, { 6, 8, 9 }, { 6, 7, 10 }, { 6, 8, 9 }, {}, {}, {}, {}, {} },
    { 6, 7, 8, 9, 10 },
    0u,

    0u,
    false,
    4u,
    -1.0,
    { { 3u, 6u } },
    {},
    {} };

TEST_CASE( "Computing gate cost", "[aqfp_network_cost]" )
{
  mockturtle::simple_cost_computer<mockturtle::aqfp_logical_network_t<int>> simp_cc( { { 3u, 3.0 } } );

  net1.computed_cost = -1.0;
  net2.computed_cost = -1.0;

  CHECK( simp_cc( net1 ) == 12.0 );
  CHECK( simp_cc( net2 ) == 18.0 );
}

TEST_CASE( "Computing AQFP cost", "[aqfp_network_cost]" )
{
  mockturtle::aqfp_cost_computer<mockturtle::aqfp_logical_network_t<int>> aqfp_cc( { { 3u, 3.0 } }, { { 3u, 3.0 } }, 1.0, 4u );

  net1.computed_cost = -1.0;
  net2.computed_cost = -1.0;

  CHECK( aqfp_cc( net1 ) == 18.0 );
  CHECK( aqfp_cc( net2 ) == 42.0 );
}
