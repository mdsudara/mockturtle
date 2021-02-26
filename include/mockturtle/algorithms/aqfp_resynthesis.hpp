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
  \file aqfp_resynthesis.hpp
  \brief Resynthesis of path balanced networks

  \author Dewmini Sudara
*/

#pragma once

#include <iostream>
#include <unordered_map>

#include "../traits.hpp"
#include "../utils/node_map.hpp"
#include "../utils/stopwatch.hpp"
#include "../views/fanout_view.hpp"
#include "../views/topo_view.hpp"

#include <fmt/format.h>

namespace mockturtle
{

/*! \brief Parameters for aqfp_resynthesis.
 *
 * The data structure `aqfp_resynthesis_params` holds configurable parameters
 * with default arguments for `node_resynthesis`.
 */
struct aqfp_resynthesis_params
{
  /*! \brief Be verbose. */
  bool verbose{ false };
  double gate_cost = 6.0;
  double buffer_cost = 2.0;
  std::map<uint32_t, double> splitters = { { 4u, 2.0 } };
};

/*! \brief Statistics of aqfp_resynthesis.
 *
 * The data structure `aqfp_resynthesis_stats` holds data collected during aqfp_resynthesis.  
 */
struct aqfp_resynthesis_stats
{
  stopwatch<>::duration time_total{ 0 };

  void report() const
  {
    std::cout << fmt::format( "[i] total time = {:>8.2f} secs\n", to_seconds( time_total ) );
  }
};

/*! \brief Results of aqfp_resynthesis.
 *
 * The data structure `aqfp_resynthesis_result` holds the resulting level assignment of nodes 
 * and the level of the critical primary output.
 */
template<typename NtkDest>
struct aqfp_resynthesis_result
{
  node_map<uint32_t, NtkDest, std::unordered_map<node<NtkDest>, uint32_t>> node_level;
  uint32_t critical_po_level = 0u;
  double gate_cost = 0.0;
  double fanout_net_cost = 0.0;
  double total_cost = 0.0;

  aqfp_resynthesis_result( NtkDest& ntk_dest ) : node_level( ntk_dest )
  {
  }

  aqfp_resynthesis_result& operator=( aqfp_resynthesis_result& other )
  {
    total_cost = other.total_cost;
    critical_po_level = other.critical_po_level;
    return *this;
  }
};

namespace detail
{

template<class NtkDest, class NtkSource, class ResynthesisFn>
class aqfp_resynthesis_impl
{
public:
  aqfp_resynthesis_impl(
      NtkDest& ntk_dest,
      NtkSource const& ntk,
      ResynthesisFn&& resynthesis_fn,
      aqfp_resynthesis_params const& ps,
      aqfp_resynthesis_result<NtkDest>& result,
      aqfp_resynthesis_stats& st )
      : ntk_dest( ntk_dest ),
        ntk( ntk ),
        resynthesis_fn( resynthesis_fn ),
        ps( ps ),
        result( result ),
        st( st )
  {
  }

  void run()
  {
    stopwatch t( st.time_total );

    node_map<signal<NtkDest>, NtkSource> node2new( ntk );

    /* map constants */
    auto c0 = ntk_dest.get_constant( false );
    node2new[ntk.get_node( ntk.get_constant( false ) )] = c0;
    result.node_level[ntk_dest.get_node( c0 )] = 0u;

    if ( ntk.get_node( ntk.get_constant( true ) ) != ntk.get_node( ntk.get_constant( false ) ) )
    {
      auto c1 = ntk_dest.get_constant( true );
      node2new[ntk.get_node( ntk.get_constant( true ) )] = c1;
      result.node_level[ntk_dest.get_node( c1 )] = 0u;
    }

    /* map primary inputs */
    ntk.foreach_pi( [&]( auto n ) {
      auto pi = ntk_dest.create_pi();
      node2new[n] = pi;
      result.node_level[ntk_dest.get_node( pi )] = 0u;

      if constexpr ( has_has_name_v<NtkSource> && has_get_name_v<NtkSource> && has_set_name_v<NtkDest> )
      {
        if ( ntk.has_name( ntk.make_signal( n ) ) )
          ntk_dest.set_name( node2new[n], ntk.get_name( ntk.make_signal( n ) ) );
      }
    } );

    ntk.foreach_ro( [&]( auto n ) {
      auto ro = ntk_dest.create_ro();
      node2new[n] = ro;
      result.node_level[ntk_dest.get_node( ro )] = 0u;

      ntk_dest._storage->latch_information[ntk_dest.get_node( node2new[n] )] = ntk._storage->latch_information[n];
      if constexpr ( has_has_name_v<NtkSource> && has_get_name_v<NtkSource> && has_set_name_v<NtkDest> )
      {
        if ( ntk.has_name( ntk.make_signal( n ) ) )
          ntk_dest.set_name( node2new[n], ntk.get_name( ntk.make_signal( n ) ) );
      }
    } );

    /* map nodes */
    fanout_view ntk_fanout{ ntk };
    depth_view ntk_depth{ ntk_fanout };
    topo_view ntk_topo{ ntk_depth };

    ntk_topo.foreach_node( [&]( auto n ) {
      if ( ntk.is_constant( n ) || ntk.is_ci( n ) )
        return;

      std::vector<signal<NtkDest>> children;
      ntk.foreach_fanin( n, [&]( auto const& f ) {
        children.push_back( ntk.is_complemented( f ) ? ntk_dest.create_not( node2new[f] ) : node2new[f] );
      } );

      bool performed_resyn = false;
      resynthesis_fn( ntk_topo, n, ntk_dest, ntk.node_function( n ), children.begin(), children.end(), result, [&]( auto const& f ) {
        node2new[n] = f;

        if constexpr ( has_has_name_v<NtkSource> && has_get_name_v<NtkSource> && has_set_name_v<NtkDest> )
        {
          if ( ntk.has_name( ntk.make_signal( n ) ) )
            ntk_dest.set_name( f, ntk.get_name( ntk.make_signal( n ) ) );
        }

        performed_resyn = true;
        return false;
      } );

      if ( !performed_resyn )
      {
        fmt::print( "[e] could not perform resynthesis for node {} in node_resynthesis\n", ntk.node_to_index( n ) );
        std::abort();
      }
    } );

    /* map primary outputs */
    ntk.foreach_po( [&]( auto const& f, auto index ) {
      (void)index;

      auto const o = ntk.is_complemented( f ) ? !node2new[f] : node2new[f];
      ntk_dest.create_po( o );

      if constexpr ( has_has_output_name_v<NtkSource> && has_get_output_name_v<NtkSource> && has_set_output_name_v<NtkDest> )
      {
        if ( ntk.has_output_name( index ) )
        {
          ntk_dest.set_output_name( index, ntk.get_output_name( index ) );
        }
      }
    } );

    ntk.foreach_ri( [&]( auto const& f, auto index ) {
      (void)index;

      auto const o = ntk.is_complemented( f ) ? ntk_dest.create_not( node2new[f] ) : node2new[f];
      ntk_dest.create_ri( o );

      if constexpr ( has_has_output_name_v<NtkSource> && has_get_output_name_v<NtkSource> && has_set_output_name_v<NtkDest> )
      {
        if ( ntk.has_output_name( index ) )
        {
          ntk_dest.set_output_name( index + ntk.num_pos(), ntk.get_output_name( index + ntk.num_pos() ) );
        }
      }
    } );
  }

  void update_cost()
  {
    fanout_view dest_fv{ ntk_dest };

    std::vector<node<NtkDest>> internal_nodes;
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
      result.gate_cost += ps.gate_cost;

      std::vector<uint32_t> rellev;

      dest_fv.foreach_fanout( n, [&]( auto fo ) {
        assert( result.node_level[fo] > result.node_level[n] );
        rellev.push_back( result.node_level[fo] - result.node_level[n] );
      } );

      uint32_t pos = 0u;
      while ( rellev.size() < dest_fv.fanout_size( n ) )
      {
        pos++;
        rellev.push_back( result.critical_po_level - result.node_level[n] );
      }

      if ( rellev.size() > 1u || ( rellev.size() == 1u && rellev[0] > 0 ) )
      {
        std::sort( rellev.begin(), rellev.end() );
        auto cst = cost_for_config( rellev );

        // if ( cst > 1000 )
        // {
        //   dest_fv.foreach_fanout( n, [&]( auto fo ) {
        //     fmt::print( "fanout {}\n", fo );
        //   } );

        //   fmt::print( "node = {}\n", n );
        //   fmt::print( "cost of rellev = {}\n", cst );
        //   fmt::print( "rellev = {}\n", fmt::join( rellev, " " ) );
        // }

        result.fanout_net_cost += cst;
      }
    }

    result.total_cost = result.gate_cost + result.fanout_net_cost;
  }

private:
  NtkDest& ntk_dest;
  NtkSource const& ntk;
  ResynthesisFn&& resynthesis_fn;
  aqfp_resynthesis_params const& ps;
  aqfp_resynthesis_result<NtkDest>& result;
  aqfp_resynthesis_stats& st;
  std::map<std::vector<uint32_t>, double> rellev_cache;

  static constexpr double IMPOSSIBLE = std::numeric_limits<double>::infinity();

  /**
   * \brief Compute the best splitter and buffer cost for a given relative level configuration 'config'.
   */
  double cost_for_config( const std::vector<uint32_t>& config )
  {
    if ( config.empty() )
      return 0.0;

    if ( config.size() == 1 )
    {
      if ( config[0] >= 1 )
      {
        return ( config[0] - 1 ) * ps.buffer_cost;
      }
      else
      {
        return IMPOSSIBLE;
      }
    }

    if ( rellev_cache.count( config ) )
    {
      return rellev_cache[config];
    }

    auto result = IMPOSSIBLE;

    for ( const auto& s : ps.splitters )
    {
      for ( auto size = 2u; size <= std::min( s.first, uint32_t( config.size() ) ); size++ )
      {
        auto sp_lev = config[config.size() - size] - 1;
        if ( sp_lev == 0 )
        {
          continue;
        }

        auto temp = s.second;

        for ( auto i = config.size() - size; i < config.size(); i++ )
        {
          temp += ( config[i] - config[config.size() - size] ) * ps.buffer_cost;
        }

        std::vector<uint32_t> new_config( config.begin(), config.begin() + ( config.size() - size ) );
        new_config.push_back( sp_lev );
        std::sort( new_config.begin(), new_config.end() );

        temp += cost_for_config( new_config );

        if ( temp < result )
        {
          result = temp;
        }
      }
    }

    return ( rellev_cache[config] = result );
  }
};

} /* namespace detail */

/*! \brief Path balanced resynthesis algorithm. */
template<class NtkDest, class NtkSource, class ResynthesisFn>
aqfp_resynthesis_result<NtkDest> aqfp_resynthesize( NtkDest& ntk_dest, NtkSource const& ntk, ResynthesisFn&& resynthesis_fn, aqfp_resynthesis_params const& ps = { false, 6.0, 2.0, { { 4u, 2.0 } } }, aqfp_resynthesis_stats* pst = nullptr )
{
  static_assert( is_network_type_v<NtkSource>, "NtkSource is not a network type" );
  static_assert( is_network_type_v<NtkDest>, "NtkDest is not a network type" );

  static_assert( has_get_node_v<NtkSource>, "NtkSource does not implement the get_node method" );
  static_assert( has_get_constant_v<NtkSource>, "NtkSource does not implement the get_constant method" );
  static_assert( has_foreach_pi_v<NtkSource>, "NtkSource does not implement the foreach_pi method" );
  static_assert( has_foreach_node_v<NtkSource>, "NtkSource does not implement the foreach_node method" );
  static_assert( has_is_constant_v<NtkSource>, "NtkSource does not implement the is_constant method" );
  static_assert( has_is_pi_v<NtkSource>, "NtkSource does not implement the is_pi method" );
  static_assert( has_is_complemented_v<NtkSource>, "NtkSource does not implement the is_complemented method" );
  static_assert( has_foreach_fanin_v<NtkSource>, "NtkSource does not implement the foreach_fanin method" );
  static_assert( has_node_function_v<NtkSource>, "NtkSource does not implement the node_function method" );
  static_assert( has_foreach_po_v<NtkSource>, "NtkSource does not implement the foreach_po method" );

  static_assert( has_get_constant_v<NtkDest>, "NtkDest does not implement the get_constant method" );
  static_assert( has_create_pi_v<NtkDest>, "NtkDest does not implement the create_pi method" );
  static_assert( has_create_not_v<NtkDest>, "NtkDest does not implement the create_not method" );
  static_assert( has_create_po_v<NtkDest>, "NtkDest does not implement the create_po method" );

  aqfp_resynthesis_result<NtkDest> result( ntk_dest );
  aqfp_resynthesis_stats st;
  resynthesis_fn.clear_state();
  detail::aqfp_resynthesis_impl<NtkDest, NtkSource, ResynthesisFn> p( ntk_dest, ntk, resynthesis_fn, ps, result, st );
  p.run();
  p.update_cost();
  if ( ps.verbose )
  {
    st.report();
  }

  if ( pst )
  {
    *pst = st;
  }

  return result;
}

} // namespace mockturtle
