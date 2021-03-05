#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <thread>
#include <vector>

#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_db_builder.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/dag_cost.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/gen_dag.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/generate_db.hpp>

#include <fmt/format.h>

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

std::unordered_map<uint32_t, uint32_t> string_to_uint_uint_map( std::string str )
{
  std::stringstream ss( str );
  std::unordered_map<uint32_t, uint32_t> res;
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

  auto num_threads = std::thread::hardware_concurrency();
  if ( num_threads == 0u )
    num_threads = 1u;
  std::cerr << fmt::format( "Will be using {} threads\n", num_threads );

  mockturtle::dag_generator_params params;
  params.allowed_num_fanins = { 3u };
  params.max_gates_of_fanin = { { 3u, 7u } };
  params.max_gates = 7u;
  params.max_levels = 7u;
  params.max_num_in = 4u;

  params.verbose = 1u;

  std::unordered_map<uint32_t, double> gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters = { { 1u, 2.0 }, { 4u, 2.0 } };

  std::string cmd( argv[1] );
  if ( cmd == "generate-dags" )
  {

    if ( argc > 2 )
    {
      if ( argc != 7 )
      {
        std::cerr << "Not enough arguments!";
        return 0;
      }

      params.allowed_num_fanins = string_to_uint_vec( std::string( argv[2] ) );
      params.max_gates_of_fanin = string_to_uint_uint_map( std::string( argv[3] ) );
      params.max_gates = std::stoul( std::string( argv[4] ) );
      params.max_levels = std::stoul( std::string( argv[5] ) );
      params.max_num_in = std::stoul( std::string( argv[6] ) );
    }

    mockturtle::generate_aqfp_dags( params, "./aqfp_dags", num_threads );
  }
  else if ( cmd == "compute-costs" )
  {
    mockturtle::compute_aqfp_dag_costs( gate_costs, splitters, "./aqfp_dags", "./aqfp_costs", num_threads );
  }
  else if ( cmd == "generate-db" )
  {
    assert( argc == 5 );

    std::string dpath( "./aqfp_dags" );
    std::string cpath( "./aqfp_costs" );
    std::string opath( "./aqfp_db" );

    mockturtle::generate_aqfp_db( gate_costs, splitters, cpath, dpath, opath, num_threads );
  }
  else if ( cmd == "db-from-scratch" )
  {
    if ( argc > 2 )
    {
      if ( argc != 7 )
      {
        std::cerr << "Not enough arguments!";
        return 0;
      }

      params.allowed_num_fanins = string_to_uint_vec( std::string( argv[2] ) );
      params.max_gates_of_fanin = string_to_uint_uint_map( std::string( argv[3] ) );
      params.max_gates = std::stoul( std::string( argv[4] ) );
      params.max_levels = std::stoul( std::string( argv[5] ) );
      params.max_num_in = std::stoul( std::string( argv[6] ) );
    }

    mockturtle::generate_aqfp_db( params, gate_costs, splitters, "./aqfp", num_threads );
  }

  else
  {
    std::cerr << fmt::format( "Invalid command {}. Must be one of the following:\n"
                              "\tgenerate-dags   -- for generating DAGs\n"
                              "\tcompute-costs   -- for costing DAGs\n"
                              "\tgenerate-db     -- for generating the AQFP database\n"
                              "\tdb-from-scratch -- for generating the AQFP database from scratch\n",
                              cmd );
  }

  return 0;
}
