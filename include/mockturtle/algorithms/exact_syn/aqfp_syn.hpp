#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>

#include <lorina/lorina.hpp>

#include "../../utils/node_map.hpp"

#include "aqfp_db.hpp"
#include "dag.hpp"

namespace mockturtle
{

/**
 * Compute the relative level configurations assuming an almost balanced splitter tree
 * for a given number of fanouts.
 */
struct splitter_counter
{
  static constexpr double IMPOSSIBLE = std::numeric_limits<double>::infinity();

  splitter_counter( const std::map<uint32_t, double>& splitters, double buffer_cost ) : splitters( splitters ), buffer_cost( buffer_cost )
  {
    // TODO: support more than one splitter type
    assert( splitters.size() == 1u );
  }

  std::tuple<uint32_t, double, std::vector<uint32_t>> balanced_splitter_tree( uint32_t num_fanouts )
  {
    if ( num_fanouts == 1u )
      return { 0u, 0.0, { 0u } };

    uint32_t sfo = splitters.begin()->first;

    uint32_t num_splitters = 1u;
    uint32_t lev = 1u;
    uint32_t t = sfo;

    while ( t < num_fanouts )
    {
      num_splitters += t;
      t *= sfo;
      lev += 1u;
    }

    // we need at lest `lev` levels, but we might not need all

    std::vector<uint32_t> levels_after( num_fanouts, lev );

    uint32_t i = 0u;
    while ( t >= num_fanouts + ( sfo - 1 ) )
    {
      num_splitters--;
      t -= ( sfo - 1 );
      levels_after[i++]--;
    }

    return { num_splitters, num_splitters * splitters.begin()->second, levels_after };
  };

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
        return ( config[0] - 1 ) * buffer_cost;
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

    for ( const auto& s : splitters )
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
          temp += ( config[i] - config[config.size() - size] ) * buffer_cost;
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

private:
  std::map<uint32_t, double> splitters;
  double buffer_cost;
  std::map<std::vector<uint32_t>, double> rellev_cache;
};

struct aqfp_syn_result
{
  double total_cost;
  double replacement_cost;
  double fanout_splitter_cost;
  double pi_balancing_cost;
  double po_balancing_cost;
  uint32_t max_level;
};

enum class aqfp_node_resyn_strategy
{
  best_cost,
  best_level
};

struct aqfp_node_resyn_param
{
  std::string db_path;
  std::map<uint32_t, double> splitters;
  double buffer_cost;
  aqfp_node_resyn_strategy strategy = aqfp_node_resyn_strategy::best_cost;
  uint32_t verbose;
};

template<typename Ntk>
struct aqfp_node_resyn
{
  aqfp_node_resyn( const aqfp_node_resyn_param& ps ) : params( ps ), db( ps.buffer_cost, ps.verbose ), splt_cntr( ps.splitters, ps.buffer_cost )
  {
    std::ifstream db_file( ps.db_path );
    assert( db_file.is_open() );
    db.load_db_from_file( db_file );
    db_file.close();
  }

  template<typename NtkSource, typename NtkDest>
  aqfp_syn_result synthesize( const NtkSource& _net, NtkDest& dest )
  {

    fanout_view net_fo_fv{ _net };
    depth_view net{ net_fo_fv };

    node_map<signal<NtkDest>, NtkSource> node2new( net );
    std::unordered_map<node<NtkDest>, uint32_t> depths;
    // node_map<uint32_t, NtkDest> depths( dest );

    std::unordered_map<uint64_t, uint32_t> node_level;
    std::map<std::tuple<uint64_t, uint64_t>, uint32_t> edge_level;
    std::vector<uint32_t> pending_po_levels;

    // get nodes of the network
    std::vector<uint64_t> nodes;
    net.foreach_node( [&]( auto n ) { nodes.push_back( n ); } );

    double replacement_cost = 0.0;
    double extra_splitter_cost = 0.0;
    uint32_t max_level = 0;

    // computing aqfp cost
    for ( auto n : nodes )
    {
      if ( net.is_pi( n ) )
      {
        const auto sig = dest.create_pi();
        node2new[n] = sig;
        node_level[n] = 0u;
        depths[dest.get_node( sig )] = 0u;
        continue;
      }

      if ( net.is_constant( n ) )
      {
        const auto sig = ( n == 0 ) ? dest.get_constant( false ) : dest.get_constant( true );
        node2new[n] = sig;
        node_level[n] = 0u;
        depths[dest.get_node( sig )] = 0u;
        continue;
      }

      std::vector<uint64_t> fanin;
      net.foreach_fanin( n, [&]( auto m ) { fanin.push_back( m ); } );

      if ( fanin.size() == 1u )
      {
        // fmt::print("FANIN ONE NODE\n");
        assert( net.node_function( n ).num_vars() == 1u );

        auto sig = node2new[net.get_node( fanin[0] )];
        if ( net.node_function( n )._bits[0] == 0x2u )
        { // buffer
          node2new[n] = net.is_complemented( fanin[0] ) ? !sig : sig;
        }
        else
        { // inverter
          node2new[n] = net.is_complemented( fanin[0] ) ? sig : !sig;
        }
        // no need to add a level
        auto const_or_pi = net.is_constant( fanin[0] ) || net.is_pi( fanin[0] );
        node_level[n] = ( const_or_pi ? node_level[fanin[0]]: edge_level[{ fanin[0], n }] + 1 );

        // depths[dest.get_node( sig )] = node_level[n];
      }
      else
      {
        while ( fanin.size() < 4u )
        {
          fanin.push_back( net.get_constant( false ) );
        }
        auto tt = kitty::extend_to( net.node_function( n ), 4u );

        std::vector<uint32_t> fanin_levels( 4u, 0u );
        std::vector<bool> fanin_is_const_or_pi( 4u, true );
        for ( auto i = 0u; i < fanin.size(); i++ )
        {
          fanin_is_const_or_pi[i] = net.is_constant( fanin[i] ) || net.is_pi( fanin[i] );
          fanin_levels[i] = ( fanin_is_const_or_pi[i] ? 0u : edge_level[{ fanin[i], n }] );
        }

        auto res = ( params.strategy == aqfp_node_resyn_strategy::best_cost ) ? db.get_best_cost_replacement( tt._bits[0], fanin_levels, fanin_is_const_or_pi ) : db.get_best_level_replacement( tt._bits[0], fanin_levels, fanin_is_const_or_pi );

        auto& repnet = std::get<1>( res ).ntk;
        auto perm = std::get<1>( res ).input_perm;
        auto invcfg = std::get<1>( res ).inverter_config;

        auto output_inv = std::get<1>( res ).output_inv;
        auto gate_levels = std::get<1>( res ).gate_levels;
        auto ind = 0u;
        std::vector<signal<NtkDest>> sigmap( repnet.nodes.size() );
        std::vector<uint32_t> levmap( repnet.nodes.size() );
        for ( auto x : repnet.input_slots )
        {
          if ( x == repnet.zero_input )
          {
            sigmap[x] = node2new[net.get_constant( false )];
            levmap[x] = 0u;
            continue;
          }
          sigmap[x] = node2new[net.get_node( fanin[perm[ind]] )];
          if ( net.is_complemented( fanin[perm[ind]] ) )
          {
            sigmap[x] = !sigmap[x];
          }
          levmap[x] = fanin_levels[perm[ind]];
          ind++;
        }

        for ( auto i = repnet.num_gates(); i > 0; i-- )
        {
          auto j = i - 1;
          auto type = invcfg[j];
          sigmap[j] = dest.create_maj(
              ( type == 0 ) ? !sigmap[repnet.nodes[j][0]] : sigmap[repnet.nodes[j][0]],
              ( type == 1 ) ? !sigmap[repnet.nodes[j][1]] : sigmap[repnet.nodes[j][1]],
              ( type == 2 ) ? !sigmap[repnet.nodes[j][2]] : sigmap[repnet.nodes[j][2]] );

          levmap[j] = 0u;
          for ( auto k = 0u; k < 3; k++ )
          {
            auto lev_dif = ( gate_levels[repnet.nodes[j][k]] > gate_levels[j] ) ? gate_levels[repnet.nodes[j][k]] - gate_levels[j] : 1u;
            levmap[j] = std::max( levmap[j], levmap[repnet.nodes[j][k]] + lev_dif );
          }
          auto current = depths[dest.get_node( sigmap[j] )];

          if ( current != 0 && current != levmap[j] )
          {
            // fmt::print( "node {} mig node {} current {} new {}\n", n, dest.get_node( sigmap[j] ), current, levmap[j] );
            // assert(false);
          }
          else
          {
            depths[dest.get_node( sigmap[j] )] = levmap[j];
          }
        }

        
        node2new[n] = output_inv ? !sigmap[0] : sigmap[0];

        if ( params.verbose > 5u )
        {
          fmt::print( "[aqfp_syn] node {} to synthesize with function {:04x}\n", n, tt._bits[0] );
          fmt::print( "[aqfp_syn] node fanin = [{}]\n", fmt::join( fanin, " " ) );
          fmt::print( "[aqfp_syn] level configuration of inputs = [{}]\n", fmt::join( fanin_levels, " " ) );
          fmt::print( "[aqfp_syn] is const or is primary input  = [{}]\n", fmt::join( fanin_is_const_or_pi, " " ) );
          fmt::print( "[aqfp_syn] can be synthesized at a cost of {} (= {} + buffer cost)\n", std::get<2>( res ), std::get<1>( res ).cost );
          fmt::print( "[aqfp_syn] level configuration of target = [{}]\n", fmt::join( std::get<1>( res ).input_levels, " " ) );
          fmt::print( "[aqfp_syn] node level after synthesis = {}\n", std::get<3>( res ) );
          fmt::print( "[aqfp_syn] network = {}\n", std::get<1>( res ).ntk.encode_as_string() );
          fmt::print( "[aqfp_syn] permutation = [{}]\n\n", fmt::join( std::get<1>( res ).input_perm, " " ) );
        }

        replacement_cost += std::get<2>( res );
        node_level[n] = std::get<3>( res );
      }

      // update fanout levels considering splitters

      auto [num_splitters, splitter_cost, splitter_levels] = splt_cntr.balanced_splitter_tree( net.fanout_size( n ) );
      (void)num_splitters;

      extra_splitter_cost += splitter_cost;

      std::vector<uint64_t> fanouts;
      net.foreach_fanout( n, [&]( auto fo ) { fanouts.push_back( fo ); } );

      std::sort( fanouts.begin(), fanouts.end(), [&]( auto f1, auto f2 ) { return ( net.depth() - net.level( f1 ) ) > ( net.depth() - net.level( f2 ) ); } );

      uint32_t foind = 0u;
      for ( auto fo : fanouts )
      {
        edge_level[{ n, fo }] = node_level[n] + splitter_levels[foind++];
      }

      for ( auto i = foind; i < net.fanout_size( n ); i++ )
      {
        pending_po_levels.push_back( node_level[n] + splitter_levels[foind++] );
      }
    }

    max_level = std::max( max_level, *std::max_element( pending_po_levels.begin(), pending_po_levels.end() ) );

    net.foreach_po( [&]( auto po ) { dest.create_po( net.is_complemented( po ) ? !node2new[net.get_node( po )] : node2new[net.get_node( po )] ); } );

    fanout_view dest_fv{ dest };

    double estimated_cost = 0.0;

    std::vector<node<mig_network>> internal_nodes;
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
      estimated_cost += 6.0;

      std::vector<uint32_t> rellev;
      dest_fv.foreach_fanout( n, [&]( auto fo ) {
        assert( depths[fo] > depths[n] );
        rellev.push_back( depths[fo] - depths[n] );
      } );

      if ( max_level > depths[n] )
      {
        while ( rellev.size() < dest_fv.fanout_size( n ) )
        {
          rellev.push_back( max_level - depths[n] );
        }
      }

      std::sort( rellev.begin(), rellev.end() );

      auto cst = splt_cntr.cost_for_config( rellev );
      estimated_cost += cst;
      if ( cst > 1000000.0 )
      {
        fmt::print( "cost inf\nrellev = {}\n", fmt::join( rellev, " " ) );
        assert( false );
      }
    }

    auto po_cost = 0.0;
    for ( auto x : pending_po_levels )
    {
      po_cost += params.buffer_cost * ( max_level - x );
    }

    auto total_cost = replacement_cost + extra_splitter_cost + po_cost;

    // fmt::print( "[aqfp-syn] estimated cost          {:8.1f}\n", estimated_cost );
    // fmt::print( "[aqfp-syn] cost of replacements    {:8.1f}\n", replacement_cost );
    // fmt::print( "[aqfp-syn] cost of extra splitters {:8.1f}\n", extra_splitter_cost );
    // fmt::print( "[aqfp-syn] cost to balance po      {:8.1f}\n", po_cost );
    // fmt::print( "[aqfp-syn] total cost              {:8.1f}\n", total_cost );

    return aqfp_syn_result{ estimated_cost, replacement_cost, extra_splitter_cost, 0.0, po_cost, max_level };
  }

private:
  aqfp_node_resyn_param params;
  aqfp_logical_network_db<aqfp_logical_network_t<int>, 4u> db;
  splitter_counter splt_cntr;
};

} // namespace mockturtle
