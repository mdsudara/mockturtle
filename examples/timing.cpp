#include <chrono>
#include <fmt/format.h>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>

#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <mockturtle/algorithms/exact_syn/dag_cost.hpp>
#include <mockturtle/algorithms/exact_syn/gen_dag.hpp>
#include <mockturtle/algorithms/exact_syn/generate_db.hpp>
#include <mockturtle/algorithms/exact_syn/sat.hpp>
 
void test_timing()
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  std::vector<uint32_t> allowed_num_fanins = { 3u };
  std::map<uint32_t, uint32_t> max_gates_of_fanin = { { 3u, 7u } };

  auto params = mockturtle::dag_generator_params();

  params.max_gates = 6u;  // allow at most 4 gates in total
  params.max_num_in = 4u; // maximum number of inputs slots (need extra one for the constant)
  params.max_levels = 3u; // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  auto trivial_cost = []( auto& net ) {(void) net; return 0.0; };
  mockturtle::dag_generator<int, decltype( trivial_cost )> gen( params, trivial_cost );

  mockturtle::aqfp_cost_computer<Ntk> aqfp_cc( { { 3u, 3.0 } }, { { 3u, 3.0 } }, 1.0, 3u );

  std::vector<uint64_t> input_tt = {
      0x0000UL,
      0xaaaaUL,
      0xccccUL,
      0xf0f0UL,
      0xff00UL,
  };
  mockturtle::dag_simulator<uint64_t> dag_sim( input_tt );

  std::vector<Ntk> all_dags;

  auto t0 = std::chrono::high_resolution_clock::now();

  auto c0 = 0u;
  while ( true )
  {
    auto x = gen.next_dag();
    if ( x == std::nullopt )
      break;

    c0++;
    if ( x.value().num_gates() == 6u )
    {
      all_dags.push_back( x.value() );
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();

  auto c1 = 0u;
  for ( auto& dag : all_dags )
  {
    aqfp_cc( dag );
    if ( ( ++c1 ) >= 1000u )
      break;
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  auto c2 = 0u;
  for ( auto& dag : all_dags )
  {
    dag_sim.all_functions_from_dag( dag );
    if ( ( ++c2 ) >= 1000u )
      break;
  }

  auto t3 = std::chrono::high_resolution_clock::now();

  std::vector<uint64_t> tts = {0x03fc, 0x19e6, 0x06b4, 0xfc26, 0x0368  };
  auto c3 = 0u;
  for (auto& dag : all_dags) {
    for (auto tt : tts) {
      mockturtle::sat_encoder se(dag, 4u, tt, 0u );
    se.solve();
    }
    if ((++c3) >= 20u) break; 
  }

  auto t4 = std::chrono::high_resolution_clock::now();

  auto d1 = std::chrono::duration_cast<std::chrono::microseconds>( t1 - t0 ).count() / 1e6;
  auto d2 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1e6;
  auto d3 = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count() / 1e6;
  auto d4 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count() / 1e6;

  fmt::print( "Generating {:5d} DAGs took                   {:5.4f} seconds\n", c0, d1 );
  fmt::print( "Computing costs for {:4d} DAGs took           {:5.4f} seconds\n", c1, d2 );
  fmt::print( "Simulating {:4d} DAGs took                    {:5.4f} seconds\n", c2, d3 );
  fmt::print( "Synthesizing {:1d} functions on {:2d} DAGs took     {:5.4f} seconds\n", tts.size(), c3, d4 );
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  test_timing();

  return 0;
}
