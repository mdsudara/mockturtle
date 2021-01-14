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
