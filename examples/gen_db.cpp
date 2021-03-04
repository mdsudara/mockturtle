#include <chrono>
#include <fmt/format.h>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>

#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_db_builder.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/dag_cost.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/gen_dag.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/generate_db.hpp>
//#include <mockturtle/algorithms/aqfp_resynthesis/sat.hpp>

void generate_dag_db( const std::vector<uint32_t>& allowed_num_fanins = { 3u }, const std::map<uint32_t, uint32_t>& max_gates_of_fanin = { { 3u, 7u } }, uint32_t max_gates = 7u, uint32_t max_num_in = 5u, uint32_t max_levels = 7u, bool count_only = false )
{
  auto params = mockturtle::dag_generator_params();

  params.max_gates = max_gates;   // allow at most 7 gates in total
  params.max_num_fanout = 1000u;  // limit the maximum fanout of a gate
  params.max_width = 1000u;       // maximum number of gates at any level
  params.max_num_in = max_num_in; // maximum number of inputs slots (need extra one for the constant)
  params.max_levels = max_levels; // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  mockturtle::generate_all_dags( params, std::cout, count_only, 1u );
}

void generate_cost_db()
{
  using Ntk = mockturtle::aqfp_dag<>;
  mockturtle::dag_aqfp_cost_all_configs<Ntk> cc( { { 3u, 6.0 } }, { { 1u, 2.0 }, { 4u, 2.0 } } );
  mockturtle::cost_all_dags( std::cin, std::cout, cc, 1u );
}

// void generate_dag_npn_db()
// {
//   mockturtle::simulate_all_dags( std::cin, std::cout, 1u );
// }

// void test_sat_based_exact_syn()
// {
//   using Ntk = mockturtle::aqfp_dag<int>;

//   uint64_t tt = 0x0001;

//   std::vector<uint32_t> allowed_num_fanins = { 3u, 5u };                        // will use only fanin 3 gates
//   std::map<uint32_t, uint32_t> max_gates_of_fanin = { { 3u, 3u }, { 5u, 0u } }; // allow at most 3 3-input gates and 3 5-input gates

//   auto params = mockturtle::dag_generator_params();

//   params.max_gates = 3u;         // allow at most 4 gates in total
//   params.max_num_fanout = 1000u; // limit the maximum fanout of a gate
//   params.max_width = 1000u;      // maximum number of gates at any level
//   params.max_num_in = 5u;        // maximum number of inputs slots (need extra one for the constant)
//   params.max_levels = 3u;        // maximum number of gate levels in a DAG

//   params.allowed_num_fanins = allowed_num_fanins;
//   params.max_gates_of_fanin = max_gates_of_fanin;

//   mockturtle::dag_generator<int, mockturtle::simple_cost_computer<Ntk>> gen( params, mockturtle::simple_cost_computer<Ntk>( { { 3u, 3.0 }, { 5u, 5.0 } } ) );
//   while ( true )
//   {
//     auto should_expand_pdag = [&]( Ntk net ) {
//       for ( auto&& t : net.last_layer_leaves )
//       {
//         net.add_leaf_node( { t } );
//       }

//       for ( auto&& t : net.other_leaves )
//       {
//         net.add_leaf_node( { t } );
//       }

//       net.last_layer_leaves.clear();
//       net.other_leaves.clear();

//       return mockturtle::is_synthesizable( net, 4u, tt, 10u );
//     };

//     auto netopt = gen.next_dag( should_expand_pdag );

//     if ( netopt == std::nullopt )
//     {
//       break;
//     }

//     auto net = netopt.value();
//     auto res = mockturtle::is_synthesizable( net, 4u, tt, 1u );

//     fmt::print( "dag {} {} synthesize 0x{:04x}\n", net.encode_as_string(), res ? "can" : "cannot", tt );
//   }
// }

void generate_aqfp_db( std::string dpath, std::string cpath, std::string opath )
{
  std::ifstream d( dpath );
  std::ifstream c( cpath );

  assert( d.is_open() && c.is_open() );

  std::string dag;
  auto dag_count = 0u;

  mockturtle::aqfp_db_builder<> db( { { 3u, 6.0 }, { 5u, 10.0 } }, { { 1u, 2.0 }, { 4u, 2.0 } } );

  while ( std::getline( d, dag ) )
  {
    dag_count++;

    // auto ok = dag_count < 10 || dag_count % 1023 == 0u; // skip many DAGs to speed-up the process
    auto ok = true; // process all DAGs

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

      mockturtle::aqfp_dag<> ntk(dag);
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

    if ( dag_count % 10000 == 0u )
    {
      std::cerr << fmt::format( "count = {}\n", dag_count );

      if ( dag_count % 100000 == 0 )
      {
        std::ofstream op( opath );
        db.save_db_to_file( op );
        op.close();
      }
    }

    if ( dag_count == 10000000 )
    {
      break;
    }
  }

  d.close();
  c.close();

  std::ofstream o( opath );
  db.save_db_to_file( o );
  o.close();
}

std::vector<uint32_t> string_to_uint_vec( std::string str )
{
  std::stringstream ss( str );
  std::vector<uint32_t> res;
  uint32_t temp;
  while ( ss >> temp )
  {
    res.push_back( temp );
  }
  return res;
}

std::map<uint32_t, uint32_t> string_to_uint_uint_map( std::string str )
{
  std::stringstream ss( str );
  std::map<uint32_t, uint32_t> res;
  uint32_t t1, t2;
  while ( ss >> t1 >> t2 )
  {
    res[t1] = t2;
  }
  return res;
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  if ( argc < 2 )
  {
    std::cerr << fmt::format( "Not enough arguments. Usage: {} cmd [opt]\n", std::string( argv[0] ) );
    return 0;
  }

  std::string cmd( argv[1] );
  if ( cmd == "gd" )
  {
    if ( argc < 3 )
    {
      generate_dag_db();
    }
    else
    {
      auto allowed_num_fanins = string_to_uint_vec( std::string( argv[2] ) );
      auto max_gates_of_fanin = string_to_uint_uint_map( std::string( argv[3] ) );
      auto max_gates = std::stoul( std::string( argv[4] ) );
      auto max_num_in = std::stoul( std::string( argv[5] ) );
      auto max_levels = std::stoul( std::string( argv[6] ) );
      generate_dag_db( allowed_num_fanins, max_gates_of_fanin, max_gates, max_num_in, max_levels, false );
    }
  }
  else if ( cmd == "cd" )
  {
    generate_cost_db();
  }
  else if ( cmd == "sd" )
  {
   // generate_dag_npn_db();
  }
  else if ( cmd == "ad" )
  {
    assert( argc == 5 );

    std::string dpath( argv[2] );
    std::string cpath( argv[3] );
    std::string opath( argv[4] );

    generate_aqfp_db( dpath, cpath, opath );
  }
  else if ( cmd == "ts" )
  {
  //  test_sat_based_exact_syn();
  }
  else
  {
    std::cerr << fmt::format( "Invalid command {}. Must be one of the following:\n"
                              "\tgd -- for generating DAGs\n"
                              "\tcd -- for costing DAGs\n"
                              "\tsd -- for simulating DAGs\n"
                              "\tad -- for generating the AQFP database\n"
                              "\tts -- for testing SAT-based exact synthesis\n",
                              cmd );
  }

  return 0;
}
