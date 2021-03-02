/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2021  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file cost_functions.hpp
  \brief Various cost functions for (optimization) algorithms

  \author Heinz Riener
  \author Mathias Soeken
  \author Dewmini Marakkalage 
*/

#pragma once

#include <cstdint>
#include <limits>

#include "../traits.hpp"

#include "../views/fanout_view.hpp"

#include "../algorithms/aqfp_resynthesis/dag_cost.hpp"

namespace mockturtle
{

template<class Ntk>
struct unit_cost
{
  uint32_t operator()( Ntk const& ntk, node<Ntk> const& node ) const
  {
    (void)ntk;
    (void)node;
    return 1u;
  }
};

template<class Ntk>
struct mc_cost
{
  uint32_t operator()( Ntk const& ntk, node<Ntk> const& node ) const
  {
    if constexpr ( has_is_xor_v<Ntk> )
    {
      if ( ntk.is_xor( node ) )
      {
        return 0u;
      }
    }

    if constexpr ( has_is_xor3_v<Ntk> )
    {
      if ( ntk.is_xor3( node ) )
      {
        return 0u;
      }
    }

    if constexpr ( has_is_nary_and_v<Ntk> )
    {
      if ( ntk.is_nary_and( node ) )
      {
        if ( ntk.fanin_size( node ) > 1u )
        {
          return ntk.fanin_size( node ) - 1u;
        }
        return 0u;
      }
    }

    if constexpr ( has_is_nary_or_v<Ntk> )
    {
      if ( ntk.is_nary_or( node ) )
      {
        if ( ntk.fanin_size( node ) > 1u )
        {
          return ntk.fanin_size( node ) - 1u;
        }
        return 0u;
      }
    }

    if constexpr ( has_is_nary_xor_v<Ntk> )
    {
      if ( ntk.is_nary_xor( node ) )
      {
        return 0u;
      }
    }

    // TODO (Does not take into account general node functions)
    return 1u;
  }
};

template<class Ntk, class NodeCostFn = unit_cost<Ntk>>
uint32_t costs( Ntk const& ntk )
{
  static_assert( is_network_type_v<Ntk>, "Ntk is not a network type" );
  static_assert( has_foreach_gate_v<Ntk>, "Ntk does not implement the foreach_gate method" );

  uint32_t total{ 0u };
  NodeCostFn cost_fn{};
  ntk.foreach_gate( [&]( auto const& n ) {
    total += cost_fn( ntk, n );
  } );
  return total;
}

struct aqfp_network_cost
{
  static constexpr double IMPOSSIBLE = std::numeric_limits<double>::infinity();

  aqfp_network_cost( const std::unordered_map<uint32_t, double>& gate_costs, const std::unordered_map<uint32_t, double>& splitters )
      : gate_costs( gate_costs ), fanout_cc( splitters ) {}

  template<typename Ntk, typename LevelMap>
  double operator()( const Ntk& ntk, const LevelMap& level_of_node, uint32_t critical_po_level )
  {
    fanout_view dest_fv{ ntk };
    auto gate_cost = 0.0;
    auto fanout_net_cost = 0.0;

    std::vector<node<Ntk>> internal_nodes;
    dest_fv.foreach_node( [&]( auto n ) {
      if ( dest_fv.is_constant( n ) || dest_fv.is_pi( n ) )
      {
        return;
      }

      if ( n > 0u && dest_fv.is_maj( n ) )
      {
        internal_nodes.push_back( n );
      }
    } );

    for ( auto n : internal_nodes )
    {
      gate_cost += gate_costs.at( ntk.fanin_size( n ) );

      std::vector<uint32_t> rellev;

      dest_fv.foreach_fanout( n, [&]( auto fo ) {
        assert( level_of_node.at( fo ) > level_of_node.at( n ) );
        rellev.push_back( level_of_node.at( fo ) - level_of_node.at( n ) );
      } );

      uint32_t pos = 0u;
      while ( rellev.size() < dest_fv.fanout_size( n ) )
      {
        pos++;
        rellev.push_back( critical_po_level - level_of_node.at( n ) );
      }

      if ( rellev.size() > 1u || ( rellev.size() == 1u && rellev[0] > 0 ) )
      {
        std::sort( rellev.begin(), rellev.end() );
        fanout_net_cost += fanout_cc( rellev );
      }
    }

    return gate_cost + fanout_net_cost;

    return 0.0;
  }

private:
  std::unordered_map<uint32_t, double> gate_costs;
  fanout_cost_computer fanout_cc;
};

} /* namespace mockturtle */