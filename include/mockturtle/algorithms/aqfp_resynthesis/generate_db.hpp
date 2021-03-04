#pragma once

#include <iostream>

#include <kitty/kitty.hpp>

#include "./dag.hpp"
#include "./dag_cost.hpp"
#include "./gen_dag.hpp"
// #include "./sat.hpp"
// #include "./simulate_dag.hpp"

namespace mockturtle
{

/**
 * \brief Generates all DAGs matching the constraints in `params`. If `count_only` is set, only the final DAG count will be printed.
 */
void generate_all_dags( const mockturtle::dag_generator_params& params, std::ostream& os, bool count_only = false, uint32_t verbose = 0u )
{
  auto gen = mockturtle::dag_generator<int, std::function<double( aqfp_dag_builder<>& )>>( params, []( aqfp_dag_builder<>& net ) { (void) net; return 0.0; } );

  auto t0 = std::chrono::high_resolution_clock::now();
  uint64_t count = 0u;
  std::vector<uint64_t> counts_inp(6);
  while ( true )
  {
    auto result_opt = gen.next_dag( []( auto& net ) {(void) net; return true; } );

    if ( result_opt == std::nullopt )
    {
      /* No more DAGs */
      if ( verbose > 0 )
      {
        std::cerr << "Finished generating DAGs" << std::endl;
      }
      break;
    }
    auto [result, cost] = result_opt.value();
    (void) cost;
    if ( !count_only )
    {
      os << fmt::format( "{}\n", result.encode_as_string() );
    }
    count++;
    counts_inp[result.input_slots.size()] ++;

    if ( verbose > 5u || ( verbose > 0 && count % 1000 == 0 ) )
    {
      auto t1 = std::chrono::high_resolution_clock::now();
      auto d1 = std::chrono::duration_cast<std::chrono::milliseconds>( t1 - t0 );
      std::cerr << fmt::format( "Number of DAGs generated {:8d}\nTime so far in seconds {:9.3f}\n", count, d1.count() / 1000.0 );
    }
  }
  auto t2 = std::chrono::high_resolution_clock::now();
  auto d2 = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t0 );
  std::cerr << fmt::format( "Number of DAGs generated {:10d}\nTime elapsed in seconds {:9.3f}\n", count, d2.count() / 1000.0 );
  std::cerr << fmt::format( "Number of DAGs of different input counts: [{}]\n", fmt::join(counts_inp, " "));
}

/**
 * \brief Compute costs for DAG considering multiple input level configurations.
 */
template<typename CostComputerT>
void cost_all_dags( std::istream& is, std::ostream& os, CostComputerT&& cc, uint32_t verbose = 0u )
{
  std::string temp;
  uint32_t count = 0u;

  while ( getline( is, temp ) )
  {
    if ( temp.length() > 0 )
    {
      if ( verbose > 5u || ( verbose > 0 && ( ++count ) % 1000u == 0u ) )
      {
        std::cerr << fmt::format( "Processing dag {} [{}]\n", ++count, temp );
      }

      mockturtle::aqfp_dag<> net(temp);

      auto costs = cc( net );

      os << "begin" << std::endl;
      for ( auto it = costs.begin(); it != costs.end(); it++ )
      {
        os << fmt::format( "{:08x} {}\n", it->first, it->second );
      }
      os << "end" << std::endl;
    }
  }
}

// /**
//  * \brief Simulates DAGs read from a file.
//  * ! Only for 4-input networks or 5-input ones with designated constant node !
//  * Reads DAG structures from the "input_stream", compute the set of
//  * npn-classes synthesizable using each DAG, and writes computed output
//  * as binary flags to "output_stream".
//  */
// void simulate_all_dags( std::istream& is, std::ostream& os, uint32_t verbose = 0u )
// {
//   std::vector<uint32_t> tt_to_id;
//   std::vector<uint64_t> id_to_npn;

//   if ( verbose > 0 )
//   {
//     std::cerr << fmt::format( "find map from truthtable to npn id...\n" );
//   }
//   detail::compute_tt_to_npn_class_mapping( 4u, tt_to_id, id_to_npn );
//   if ( verbose > 0 )
//   {
//     std::cerr << fmt::format( "found map from truthtable to npn id...\n" );
//   }

//   std::vector<uint64_t> input_tt = {
//       0x0000UL,
//       0xaaaaUL,
//       0xccccUL,
//       0xf0f0UL,
//       0xff00UL,
//   };

//   dag_simulator<uint64_t> dag_sim( input_tt );

//   std::string temp;
//   uint32_t count = 0u;
//   while ( getline( is, temp ) )
//   {
//     if ( temp.length() > 0 )
//     {
//       if ( verbose > 5u || ( verbose > 0 && ( ++count ) % 1000 == 0 ) )
//       {
//         std::cerr << fmt::format( "processing dag {} [{}]\n", count, temp );
//       }

//       mockturtle::aqfp_logical_network_t<int> net;
//       net.decode_dag( temp );

//       if ( net.zero_input == 0 && net.input_slots.size() > 4u )
//       {
//         continue;
//       }

//       uint64_t npn_flags[4] = { 0ul, 0ul, 0ul, 0ul };
//       auto funcs = dag_sim.all_functions_from_dag( net );
//       for ( auto f : funcs )
//       {
//         auto id = tt_to_id[f];
//         npn_flags[id / 64ul] |= ( 1ul << ( id % 64ul ) );
//       }

//       os << fmt::format( "{:16x} {:16x} {:16x} {:16x}\n", npn_flags[0], npn_flags[1], npn_flags[2], npn_flags[3] );
//     }
//   }
// }

// template<typename Ntk, typename TruthTableT>
// bool is_synthesizable( const Ntk& ntk, const int32_t num_inputs, const TruthTableT& tt, uint32_t verbose )
// {
//   sat_encoder se( ntk, num_inputs, tt, verbose );
//   return se.solve();
// }

} // namespace mockturtle
