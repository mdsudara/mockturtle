#include <chrono>
#include <fmt/format.h>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>
#include <fstream>

#include <mockturtle/algorithms/exact_syn/aqfp_db.hpp>
#include <mockturtle/algorithms/exact_syn/dag_cost.hpp>
#include <mockturtle/algorithms/exact_syn/gen_dag.hpp>
#include <mockturtle/algorithms/exact_syn/generate_db.hpp>
#include <mockturtle/algorithms/exact_syn/sat.hpp>

void generate_dag_db()
{
  std::vector<uint32_t> allowed_num_fanins = { 3u };                // will use only fanin 3 gates
  std::map<uint32_t, uint32_t> max_gates_of_fanin = { { 3u, 7u } }; // allow at most 7 3-input gates and 0 5-input gates

  auto params = mockturtle::dag_generator_params();

  params.max_gates = 7u;         // allow at most 7 gates in total
  params.max_num_fanout = 1000u; // limit the maximum fanout of a gate
  params.max_width = 1000u;      // maximum number of gates at any level
  params.max_num_in = 5u;        // maximum number of inputs slots (need extra one for the constant)
  params.max_levels = 7u;        // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  mockturtle::generate_all_dags( params, std::cout, 1u );
}

void generate_cost_db()
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;
  mockturtle::aqfp_cost_computer<Ntk> cc( { { 3u, 3.0 } }, { { 3u, 3.0 } }, 1.0, 4u );
  mockturtle::cost_all_dags( std::cin, std::cout, cc, 1u );
}

void generate_dag_npn_db()
{
  mockturtle::simulate_all_dags( std::cin, std::cout, 1u );
}

void test_sat_based_exact_syn()
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  uint64_t tt = 0x0001;

  std::vector<uint32_t> allowed_num_fanins = { 3u, 5u };                        // will use only fanin 3 gates
  std::map<uint32_t, uint32_t> max_gates_of_fanin = { { 3u, 3u }, { 5u, 0u } }; // allow at most 3 3-input gates and 3 5-input gates

  auto params = mockturtle::dag_generator_params();

  params.max_gates = 3u;         // allow at most 4 gates in total
  params.max_num_fanout = 1000u; // limit the maximum fanout of a gate
  params.max_width = 1000u;      // maximum number of gates at any level
  params.max_num_in = 5u;        // maximum number of inputs slots (need extra one for the constant)
  params.max_levels = 3u;        // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  mockturtle::dag_generator<int, mockturtle::simple_cost_computer<Ntk>> gen( params, mockturtle::simple_cost_computer<Ntk>( { { 3u, 3.0 }, { 5u, 5.0 } } ) );
  while ( true )
  {
    auto should_expand_pdag = [&]( Ntk net ) {
      for ( auto&& t : net.last_layer_leaves )
      {
        net.add_leaf_node( { t } );
      }

      for ( auto&& t : net.other_leaves )
      {
        net.add_leaf_node( { t } );
      }

      net.last_layer_leaves.clear();
      net.other_leaves.clear();

      return mockturtle::is_synthesizable( net, 4u, tt, 10u );
    };

    auto netopt = gen.next_dag( should_expand_pdag );

    if ( netopt == std::nullopt )
    {
      break;
    }

    auto net = netopt.value();
    auto res = mockturtle::is_synthesizable( net, 4u, tt, 1u );

    fmt::print( "dag {} {} synthesize 0x{:04x}\n", net.encode_as_string(), res ? "can" : "cannot", tt );
  }
}

void generate_aqfp_db( std::string dpath, std::string cpath, std::string opath )
{
  std::ifstream d( dpath );
  std::ifstream c( cpath );

  assert (d.is_open() && c.is_open());

  std::string dag;
  auto dag_count = 0u;

  mockturtle::aqfp_logical_network_db<mockturtle::aqfp_logical_network_t<int>, 4u> db( 1u );

  while ( std::getline( d, dag ) )
  {
    dag_count++;
    auto ok = dag_count < 10 || dag_count % 1023 == 0u;

    if ( ok )
    {
      std::string token;
      c >> token;
      std::unordered_map<uint64_t, double> configs;
      while ( true )
      {
        double cst;
        c >> token;
        if ( token == "end" )
        {
          break;
        }
        c >> cst;
        configs[std::stoul( token, 0, 16 )] = cst;
      }

      mockturtle::aqfp_logical_network_t<int> ntk;
      ntk.decode_dag( dag );
      if ( ntk.input_slots.size() < 5u || ( ntk.input_slots.size() == 5u && ntk.zero_input != 0 ) )
      {
        db.update( ntk, configs );
      }
    }
    else
    {
      std::string token;
      c >> token;
      while ( true )
      {
        c >> token;
        if ( token == "end" )
        {
          break;
        }
        c >> token;
      }
    }

    if ( dag_count % 100000 == 0u )
    {
      std::cerr << fmt::format( "count = {}\n", dag_count );
    }

    if ( dag_count == 10000000 )
    {
      break;
    }
  }
  d.close();
  c.close();

  std::ofstream o( opath);
  db.save_db_to_file( o );
  o.close();
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  if ( argc < 2 )
  {
    std::cerr << fmt::format( "Not enough arguments.\n" );
    return 0;
  }

  std::string cmd( argv[1] );
  if ( cmd == "gd" )
  {
    generate_dag_db();
  }
  else if ( cmd == "cd" )
  {
    generate_cost_db();
  }
  else if ( cmd == "sd" )
  {
    generate_dag_npn_db();
  }
  else if ( cmd == "ts" )
  {
    test_sat_based_exact_syn();
  }
  else if ( cmd == "ad" )
  {
    assert(argc == 5);

    std::string dpath(argv[2]);
    std::string cpath(argv[3]);
    std::string opath(argv[4]);

    generate_aqfp_db(dpath, cpath, opath);
  }
  else
  {
    std::cerr << fmt::format( "Invalid command {}. Must be one of the following:\n"
                              "\tgd -- for generating DAGs\n"
                              "\tcd -- for costing DAGs\n"
                              "\tsd -- for simulating DAGs\n"
                              "\tts -- for testing SAT-based exact synthesis\n",
                              cmd );
  }

  return 0;
}
