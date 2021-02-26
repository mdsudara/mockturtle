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

  splitter_counter( const std::map<uint32_t, double>& splitters ) : splitters( splitters )
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

private:
  std::map<uint32_t, double> splitters;
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

template<typename SrcNodeT, typename DestNodeT>
struct aqfp_resyn
{
  aqfp_resyn( const aqfp_node_resyn_param& ps ) : params( ps ), db( ps.buffer_cost, ps.verbose ), splt_cntr( ps.splitters )
  {
    std::ifstream db_file( ps.db_path );
    assert( db_file.is_open() );
    db.load_db_from_file( db_file );
    db_file.close();
  }

  template<typename NtkSource, typename NtkDest, typename TruthTable, typename LeavesIterator, typename Result, typename Callback>
  void operator()(
      const NtkSource& ntk_src,
      node<NtkSource> n,
      NtkDest& ntk_dest,
      const TruthTable& f,
      LeavesIterator leaves_begin,
      LeavesIterator leaves_end,
      Result& result,
      Callback&& callback )
  {

    static_assert( has_foreach_fanout_v<NtkSource>, "NtkSource does not implemente the foreach_fanout method");

    std::vector<signal<NtkDest>> fanin;

    for ( auto it = leaves_begin; it != leaves_end; it++ )
    {
      fanin.push_back( *it );
    }

    auto new_n = ntk_dest.get_constant( false );
    auto new_n_const_or_ci = false;
    auto new_n_base_level = 0u;

    if ( fanin.size() == 1u )
    {
      /* Node n is a buffer, inverter, or essentially a constant. */
      assert( f.num_vars() == 1u );

      switch ( f._bits[0] )
      {
      case 0x0u:
        new_n = ntk_dest.get_constant( false );
        new_n_const_or_ci = true;
        break;
      case 0x1u:
        new_n = !fanin[0];
        break;
      case 0x2u:
        new_n = fanin[0];
        break;
      case 0x3u:
        new_n = ntk_dest.get_constant( true );
        new_n_const_or_ci = true;
        break;
      default:
        assert( false );
      }

      new_n_const_or_ci = ntk_dest.is_constant( ntk_dest.get_node( new_n ) ) || ntk_dest.is_ci( ntk_dest.get_node( new_n ) );
      new_n_base_level = new_n_const_or_ci ? 0u : edge_level[{ ntk_dest.get_node(fanin[0]), n }] + 1;
    }
    else
    {
      // should not have more than 4 fanin nodes
      assert( fanin.size() <= 4u );

      // if less than 4 fanins, add dummy inputs
      while ( fanin.size() < 4u )
      {
        fanin.push_back( ntk_dest.get_constant( false ) );
      }
      auto tt = kitty::extend_to( f, 4u );

      std::vector<uint32_t> fanin_levels( 4u, 0u );
      std::vector<bool> fanin_const_or_ci( 4u, true );
      for ( auto i = 0u; i < fanin.size(); i++ )
      {
        fanin_const_or_ci[i] = ntk_dest.is_constant( ntk_dest.get_node( fanin[i] ) ) || ntk_dest.is_ci( ntk_dest.get_node( fanin[i] ) );
        fanin_levels[i] = ( fanin_const_or_ci[i] ? 0u : edge_level[{ ntk_dest.get_node(fanin[i]), n }] );
      }

      auto res = ( params.strategy == aqfp_node_resyn_strategy::best_cost )
                     ? db.get_best_cost_replacement( tt._bits[0], fanin_levels, fanin_const_or_ci )
                     : db.get_best_level_replacement( tt._bits[0], fanin_levels, fanin_const_or_ci );

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
          sigmap[x] = ntk_dest.get_constant( false );
          levmap[x] = 0u;
          continue;
        }

        sigmap[x] = fanin[perm[ind]];
        levmap[x] = fanin_levels[perm[ind]];
        ind++;
      }

      for ( auto i = repnet.num_gates(); i > 0; i-- )
      {
        auto j = i - 1;
        auto type = invcfg[j];
        sigmap[j] = ntk_dest.create_maj(
            ( type == 0 ) ? !sigmap[repnet.nodes[j][0]] : sigmap[repnet.nodes[j][0]],
            ( type == 1 ) ? !sigmap[repnet.nodes[j][1]] : sigmap[repnet.nodes[j][1]],
            ( type == 2 ) ? !sigmap[repnet.nodes[j][2]] : sigmap[repnet.nodes[j][2]] );

        levmap[j] = 0u;
        for ( auto k = 0u; k < repnet.nodes[j].size(); k++ )
        {
          auto lev_dif = ( gate_levels[repnet.nodes[j][k]] > gate_levels[j] ) ? gate_levels[repnet.nodes[j][k]] - gate_levels[j] : 1u;
          levmap[j] = std::max( levmap[j], levmap[repnet.nodes[j][k]] + lev_dif );
        }

        // if we have already set the level, do not increase it
        // auto current = result.node_level[ntk_dest.get_node( sigmap[j] )];
        if ( !result.node_level.has( sigmap[j] ) )
        {
          result.node_level[sigmap[j]] = levmap[j];
        }
      }

      new_n = output_inv ? !sigmap[0] : sigmap[0];
      new_n_const_or_ci = ntk_dest.is_constant( ntk_dest.get_node(new_n) ) || ntk_dest.is_ci( ntk_dest.get_node(new_n) );
      new_n_base_level = new_n_const_or_ci ? 0u : std::get<3>( res );
    }

    callback( new_n );

    // update fanout levels considering splitters
    auto [num_splitters, splitter_cost, splitter_levels] = splt_cntr.balanced_splitter_tree( ntk_src.fanout_size( n ) );
    (void)num_splitters;
    (void)splitter_cost;

    std::vector<SrcNodeT> fanouts;
    ntk_src.foreach_fanout( n, [&]( auto fo ) { fanouts.push_back( fo ); } );

    std::sort( fanouts.begin(), fanouts.end(), [&]( auto f1, auto f2 ) { return ( ntk_src.depth() - ntk_src.level( f1 ) ) > ( ntk_src.depth() - ntk_src.level( f2 ) ); } );

    uint32_t foind = 0u;
    for ( auto fo : fanouts )
    {
      if ( !new_n_const_or_ci )
      {
        edge_level[{ ntk_dest.get_node(new_n), fo }] = new_n_base_level + splitter_levels[foind++];
      }
    }

    for ( auto i = foind; i < ntk_src.fanout_size( n ); i++ )
    {
      if ( new_n_const_or_ci )
      {
        result.critical_po_level = std::max( result.critical_po_level, 0u );
      }
      else
      {
        result.critical_po_level = std::max( result.critical_po_level, new_n_base_level + splitter_levels[foind++] );
      }
    }
  }

  void clear_state() {
    edge_level.clear();
  }

private:
  aqfp_node_resyn_param params;
  aqfp_logical_network_db<aqfp_logical_network_t<int>, 4u> db;
  splitter_counter splt_cntr;

  std::map<std::tuple<DestNodeT, SrcNodeT>, uint32_t> edge_level;
};

} // namespace mockturtle
