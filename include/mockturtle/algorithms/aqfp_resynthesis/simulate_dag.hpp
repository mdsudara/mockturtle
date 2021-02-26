#pragma once

#include <unordered_set>
#include <vector>

#include <fmt/format.h>

#include "./gen_dag.hpp"

namespace mockturtle
{

namespace detail
{

/**
 * Compute the majority function on truth-tables.
 */
template<typename T>
inline T maj( const T& a, const T& b, const T& c )
{
  return ( a & b ) | ( c & ( a | b ) );
}

/**
 *  For each num_var-input function, compute its npn-class, assign each npn-class a unique id from 0 to (num-npn-classe - 1).
 *  Save the corresponding mappings in tt_to_id and id_to_npn.
 */
void compute_tt_to_npn_class_mapping( uint32_t num_vars, std::vector<uint32_t>& tt_to_id, std::vector<uint64_t>& id_to_npn )
{
  std::unordered_map<uint64_t, uint32_t> npn_to_id;

  tt_to_id = std::vector<uint32_t>( ( 1ul << ( 1ul << num_vars ) ), 0u );

  kitty::dynamic_truth_table dtt( num_vars );
  do
  {
    auto npn = kitty::exact_npn_canonization( dtt );
    auto npn_dtt = std::get<0>( npn );

    if ( !npn_to_id.count( npn_dtt._bits[0] ) )
    {
      assert( npn_dtt == dtt );
      npn_to_id[npn_dtt._bits[0]] = npn_to_id.size();
    }
    tt_to_id[dtt._bits[0]] = npn_to_id[npn_dtt._bits[0]];

    kitty::next_inplace( dtt );
  } while ( !kitty::is_const0( dtt ) );

  id_to_npn = std::vector<uint64_t>( npn_to_id.size(), 0ul );
  for ( auto it = npn_to_id.begin(); it != npn_to_id.end(); it++ )
  {
    id_to_npn[it->second] = it->first;
  }
}

} // namespace detail

template<typename TruthTableT>
class dag_simulator
{
public:
  dag_simulator( const std::vector<TruthTableT>& input_tt ) : input_tt( input_tt )
  {
  }

  std::unordered_set<TruthTableT> all_functions_from_dag( const mockturtle::aqfp_logical_network_t<int>& net )
  {
    assert( !net.is_partial_dag );

    uint32_t num_inputs = net.input_slots.size();
    if ( net.zero_input != 0 )
    {
      num_inputs--;
    }

    assert( num_inputs <= 6 );

    std::unordered_set<TruthTableT> res;

    std::vector<TruthTableT> tt( net.nodes.size(), input_tt[0] );
    auto input_ind = 1u;
    for ( auto i : net.input_slots )
    {
      if ( i != (int)net.zero_input )
      {
        tt[i] = input_tt[input_ind++];
      }
    }

    for ( auto inv_config = 0ul; inv_config < ( 1ul << ( 2 * net.num_gates() ) ); inv_config++ )
    {
      for ( auto i = net.num_gates(); i > 0; i-- )
      {
        auto ith_gate_config = ( inv_config >> ( 2 * ( i - 1 ) ) ) & 3ul;
        tt[i - 1] = detail::maj(
            ith_gate_config == 0ul ? ~tt[net.nodes[i - 1][0]] : tt[net.nodes[i - 1][0]],
            ith_gate_config == 1ul ? ~tt[net.nodes[i - 1][1]] : tt[net.nodes[i - 1][1]],
            ith_gate_config == 2ul ? ~tt[net.nodes[i - 1][2]] : tt[net.nodes[i - 1][2]] );
      }
      res.insert( tt[0] );
    }

    return res;
  }

private:
  std::vector<TruthTableT> input_tt;
};

} // namespace mockturtle
