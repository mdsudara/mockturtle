#include <catch.hpp>

#include <set>

#include <fmt/format.h>

#include <mockturtle/algorithms/exact_syn/gen_dag_util.hpp>
#include <mockturtle/algorithms/exact_syn/gen_dag.hpp>

using namespace mockturtle;

TEST_CASE( "Partition generation", "[aqfp_network_gen]" )
{
  using part = std::multiset<int>;
  using partition = std::multiset<part>;
  using partition_set = std::set<partition>;

  detail::partition_generator<int> partition_gen;

  /* General partitioning with distinct elements and unit max counts with unlimited part count and unlimited part size */
  const auto t1a = partition_gen( { 0, 1, 2, 3 }, { 1u, 1u, 1u, 1u } );
  const auto t1b = partition_gen( { 0, 1, 2, 3 }, { 1u, 1u, 1u, 1u }, 0u, 0u );
  const partition_set s1 = {
	  { { 0 }, { 1 }, { 2 }, { 3 } },
	  { { 0 }, { 1 }, { 2, 3 } },
	  { { 0 }, { 1, 2 }, { 3 } },
	  { { 0 }, { 1, 2, 3 } },
	  { { 0 }, { 1, 3 }, { 2 } },
	  { { 0, 1 }, { 2 }, { 3 } },
	  { { 0, 1 }, { 2, 3 } },
	  { { 0, 1, 2 }, { 3 } },
	  { { 0, 1, 2, 3 } },
	  { { 0, 1, 3 }, { 2 } },
	  { { 0, 2 }, { 1 }, { 3 } },
	  { { 0, 2 }, { 1, 3 } },
	  { { 0, 2, 3 }, { 1 } },
	  { { 0, 3 }, { 1 }, { 2 } },
	  { { 0, 3 }, { 1, 2 } } };
  CHECK( t1a == t1b );
  CHECK( t1a == s1 );

  /* Duplicate element and unit max counts */
  const auto t2 = partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 0u );
  const partition_set s2 = {
	  { { 0 }, { 1 }, { 1 }, { 3 } },
	  { { 0 }, { 1 }, { 1, 3 } },
	  { { 0, 1 }, { 1 }, { 3 } },
	  { { 0, 1 }, { 1, 3 } },
	  { { 0, 1, 3 }, { 1 } },
	  { { 0, 3 }, { 1 }, { 1 } } };
  CHECK( t2 == s2 );

  /* Duplicate element and unit max counts with limited part count */
  const auto t3 = partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 2u );
  const partition_set s3 = {
	  { { 0, 1 }, { 1, 3 } },
	  { { 0, 1, 3 }, { 1 } } };
  CHECK( t3 == s3 );

  /* Duplicate element and unit max counts with limited part size */
  const auto t4 = partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 0u, 2u );
  const partition_set s4 = {
	  { { 0 }, { 1 }, { 1 }, { 3 } },
	  { { 0 }, { 1 }, { 1, 3 } },
	  { { 0, 1 }, { 1 }, { 3 } },
	  { { 0, 1 }, { 1, 3 } },
	  { { 0, 3 }, { 1 }, { 1 } } };
  CHECK( t4 == s4 );

  /* Duplicate element and unit max counts with limited part count and limited part size */
  const auto t5 = partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 2u, 2u );
  const partition_set s5 = {
	  { { 0, 1 }, { 1, 3 } } };
  CHECK( t5 == s5 );

  /* Duplicate element and max count of 2 for that element */
  const auto t6 = partition_gen( { 0, 1, 1, 3 }, { 1u, 2u, 1u, 1u } );
  const partition_set s6 = {
	  { { 0 }, { 1 }, { 1 }, { 3 } },
	  { { 0 }, { 1 }, { 1, 3 } },
	  { { 0 }, { 1, 1 }, { 3 } },
	  { { 0 }, { 1, 1, 3 } },
	  { { 0, 1 }, { 1 }, { 3 } },
	  { { 0, 1 }, { 1, 3 } },
	  { { 0, 1, 1 }, { 3 } },
	  { { 0, 1, 1, 3 } },
	  { { 0, 1, 3 }, { 1 } },
	  { { 0, 3 }, { 1 }, { 1 } },
	  { { 0, 3 }, { 1, 1 } } };
  CHECK( t6 == s6 );

  /* Duplicate element and max count of 2 for that element with limited part size */
  const auto t7 = partition_gen( { 0, 1, 1, 3 }, { 1u, 2u, 1u, 1u }, 0u, 2u );
  const partition_set s7 = {
	  { { 0 }, { 1 }, { 1 }, { 3 } },
	  { { 0 }, { 1 }, { 1, 3 } },
	  { { 0 }, { 1, 1 }, { 3 } },
	  { { 0, 1 }, { 1 }, { 3 } },
	  { { 0, 1 }, { 1, 3 } },
	  { { 0, 3 }, { 1 }, { 1 } },
	  { { 0, 3 }, { 1, 1 } } };
  CHECK( t7 == s7 );
}

TEST_CASE( "Partition extension", "[aqfp_network_gen]" )
{
  using part = std::multiset<int>;
  using partition = std::multiset<part>;
  using partition_set = std::set<partition>;

  detail::partition_extender<int> partition_ext;

  auto t1 = partition_ext( { 0, 1, 2, 3, 4 }, {}, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s1 = {};
  CHECK( t1 == s1 );

  auto t2 = partition_ext( {}, { { 1, 2 }, { 3, 4 } }, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s2 = {
	  { { 1, 2 }, { 3, 4 } } };
  CHECK( t2 == s2 );

  auto t3 = partition_ext( { 0, 1, 2, 3, 4 }, { {} }, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s3 = {
	  { { 0, 1, 2, 3, 4 } } };
  CHECK( t3 == s3 );

  auto t4 = partition_ext( { 4 }, { { 1 }, { 2 }, { 3 }, { 4 } }, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s4 = {
	  { { 1 }, { 2 }, { 3, 4 }, { 4 } },
	  { { 1 }, { 2, 4 }, { 3 }, { 4 } },
	  { { 1, 4 }, { 2 }, { 3 }, { 4 } } };
  CHECK( t4 == s4 );

  auto t5 = partition_ext( { 4 }, { { 1 }, { 2 }, { 3 }, { 4 } }, { 1u, 1u, 1u, 1u, 2u } );
  partition_set s5 = {
	  { { 1 }, { 2 }, { 3 }, { 4, 4 } },
	  { { 1 }, { 2 }, { 3, 4 }, { 4 } },
	  { { 1 }, { 2, 4 }, { 3 }, { 4 } },
	  { { 1, 4 }, { 2 }, { 3 }, { 4 } } };
  CHECK( t5 == s5 );

  auto t6 = partition_ext( { 3, 4 }, { { 0 }, { 1 }, { 2 } }, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s6 = {
	  { { 0 }, { 1 }, { 2, 3, 4 } },
	  { { 0 }, { 1, 3 }, { 2, 4 } },
	  { { 0 }, { 1, 3, 4 }, { 2 } },
	  { { 0 }, { 1, 4 }, { 2, 3 } },
	  { { 0, 3 }, { 1 }, { 2, 4 } },
	  { { 0, 3 }, { 1, 4 }, { 2 } },
	  { { 0, 3, 4 }, { 1 }, { 2 } },
	  { { 0, 4 }, { 1 }, { 2, 3 } },
	  { { 0, 4 }, { 1, 3 }, { 2 } } };
  CHECK( t6 == s6 );

  auto t7 = partition_ext( { 0, 1, 2 }, { { 3 }, { 4 } }, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s7 = {
	  { { 0, 1, 2, 3 }, { 4 } },
	  { { 0, 1, 2, 4 }, { 3 } },
	  { { 0, 1, 3 }, { 2, 4 } },
	  { { 0, 1, 4 }, { 2, 3 } },
	  { { 0, 2, 3 }, { 1, 4 } },
	  { { 0, 2, 4 }, { 1, 3 } },
	  { { 0, 3 }, { 1, 2, 4 } },
	  { { 0, 4 }, { 1, 2, 3 } } };
  CHECK( t7 == s7 );

  auto t8 = partition_ext( { 3, 4, 4, 4 }, { { 0 }, { 1 }, { 2 } }, { 1u, 1u, 1u, 1u, 1u } );
  partition_set s8 = {
	  { { 0, 3, 4 }, { 1, 4 }, { 2, 4 } },
	  { { 0, 4 }, { 1, 3, 4 }, { 2, 4 } },
	  { { 0, 4 }, { 1, 4 }, { 2, 3, 4 } } };
  CHECK( t8 == s8 );
}

TEST_CASE( "Sublist generation", "[aqfp_network_gen]" )
{
  //TODO
}

TEST_CASE( "Layer extension from partial DAGs", "[aqfp_network_gen]" )
{
  //TODO
}

TEST_CASE( "DAG generation from partial DAGs", "[aqfp_network_gen]" )
{
  //TODO
}

TEST_CASE( "DAG generation with cost", "[aqfp_network_gen]" )
{
  //TODO
}

