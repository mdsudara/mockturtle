#pragma once

#include <chrono>
#include <fmt/format.h>
#include <kitty/kitty.hpp>
#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using truth_table_t = uint64_t;





uint64_t input_tt[] = {
    0x0000000000000000UL,
    0xaaaaaaaaaaaaaaaaUL,
    0xccccccccccccccccUL,
    0xf0f0f0f0f0f0f0f0UL,
    0xff00ff00ff00ff00UL,
    0xffff0000ffff0000UL,
    0xffffffff00000000UL,
};

inline uint64_t maj( uint64_t a, uint64_t b, uint64_t c )
{
  return ( a & b ) | ( c & ( a | b ) );
}

std::unordered_set<truth_table_t> all_functions_from_dag( const mockturtle::aqfp_logical_network& net )
{
  assert( net.is_dag );

  size_t num_inputs = net.input_slots.size();
  if ( net.zero_input != 0 )
    num_inputs--;

  assert( num_inputs <= 6 );

  std::unordered_set<uint64_t> res;

  std::vector<truth_table_t> tt( net.nodes.size(), input_tt[0] );
  auto input_ind = 1u;
  for ( auto i : net.input_slots )
  {
    if ( i != net.zero_input )
    {
      tt[i] = input_tt[input_ind++];
    }
  }

  for ( auto inv_config = 0ul; inv_config < ( 1ul << ( 2 * net.num_gates() ) ); inv_config++ )
  {
    for ( auto i = net.num_gates(); i > 0; i-- )
    {
      auto ith_gate_config = ( inv_config >> ( 2 * ( i - 1 ) ) ) & 3ul;
      tt[i - 1] = maj(
          ith_gate_config == 0ul ? ~tt[net.nodes[i - 1][0]] : tt[net.nodes[i - 1][0]],
          ith_gate_config == 1ul ? ~tt[net.nodes[i - 1][1]] : tt[net.nodes[i - 1][1]],
          ith_gate_config == 2ul ? ~tt[net.nodes[i - 1][2]] : tt[net.nodes[i - 1][2]] );
    }
    // fmt::print("inv config {:16x} tts = {}\n", inv_config, fmt::join(tt, " "));
    res.insert( tt[0] & 0xffff );
  }

  return res;
}

std::vector<uint64_t> compute_4_input_tt_to_npn_class_mapping()
{
  std::unordered_map<uint64_t, size_t> npn_to_id;
  std::vector<uint64_t> res( ( 1ul << ( 1ul << 4 ) ) );

  kitty::dynamic_truth_table dtt( 4u );
  do
  {
    auto npn = kitty::exact_npn_canonization( dtt );
    auto npn_dtt = std::get<0>( npn );

    if ( !npn_to_id.count( npn_dtt._bits[0] ) )
    {
      assert( npn_dtt == dtt );
      npn_to_id[npn_dtt._bits[0]] = npn_to_id.size();
    }
    res[dtt._bits[0]] = npn_to_id[npn_dtt._bits[0]];

    kitty::next_inplace( dtt );
  } while ( !kitty::is_const0( dtt ) );

  assert( npn_to_id.size() == 222ul );

  return res;
}

void simulate_all_dags()
{
  std::string temp;
  size_t count = 0u;

  fmt::print("find map from truthtable to npn id...\n");
  auto tt_to_npn_id_map = compute_4_input_tt_to_npn_class_mapping();
  fmt::print("found map from truthtable to npn id...\n");

  while ( getline( std::cin, temp ) )
  {
    if ( temp.length() > 0 )
    {
      fmt::print( "processing dag {} [{}]\n", ++count, temp );
      mockturtle::aqfp_logical_network_t net;
      net.decode_string( temp );

      uint64_t npn_flags[4] = {0ul};
      auto funcs = all_functions_from_dag(net);
      for (auto f : funcs) {
//        fmt::print("dag can synthesize tt {}\n", f);
        auto id = tt_to_npn_id_map[f];
 //       fmt::print("id for tt {}\n", id);
        npn_flags[id / 64ul] |= (1ul << (id % 64ul));
      }

      fmt::print("{:16x} {:16x} {:16x} {:16x}\n", npn_flags[0], npn_flags[1], npn_flags[2], npn_flags[3]);
    }
  }
}