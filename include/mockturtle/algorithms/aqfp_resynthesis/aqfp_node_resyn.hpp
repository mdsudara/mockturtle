#pragma once

#include "../../traits.hpp"

#include "./aqfp_db.hpp"

#include <limits>
#include <map>
#include <unordered_map>
#include <type_traits>
#include <vector>

namespace mockturtle
{

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
struct aqfp_node_resyn
{
  aqfp_node_resyn( const aqfp_node_resyn_param& ps ) : params( ps ), db( ps.buffer_cost, ps.verbose )
  {
    std::ifstream db_file( ps.db_path );
    assert( db_file.is_open() );
    db.load_db_from_file( db_file );
    db_file.close();
  }

  template<typename NtkDest, typename TruthTable, typename LeavesIterator, typename LevelUpdateCallback, typename ResynPerformedCallback>
  void operator()(
      NtkDest& ntk_dest,
      const TruthTable& f,
      LeavesIterator leaves_begin,
      LeavesIterator leaves_end,
      LevelUpdateCallback&& level_update_callback,
      ResynPerformedCallback&& resyn_performed_callback )
  {
    static_assert( std::is_invocable_v<LevelUpdateCallback, node<NtkDest>, uint32_t>, "LevelUpdateCallback must be callable with arguments of types (node<NtkDest>, level)" );
    static_assert( std::is_invocable_v<ResynPerformedCallback, signal<NtkDest>>, "ResynPerformedCallback must be callable with an argument of type signal<NtkDest>" );

    std::vector<signal<NtkDest>> leaves;
    std::vector<uint32_t> leaf_levels;
    std::vector<bool> leaf_no_splitters;
    for ( auto it = leaves_begin; it != leaves_end; it++ )
    {
      auto leaf = std::get<0>(*it);
      auto leaf_level = std::get<1>(*it);

      leaves.push_back( leaf );
      leaf_levels.push_back( leaf_level );
      leaf_no_splitters.push_back( ntk_dest.is_constant( ntk_dest.get_node( leaf ) ) || ntk_dest.is_ci( ntk_dest.get_node( leaf ) ) );
    }

    auto new_n = ntk_dest.get_constant( false );
    auto new_n_const_or_ci = false;
    auto new_n_base_level = 0u;

    if ( leaves.size() == 1u )
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
        new_n = !leaves[0];
        break;
      case 0x2u:
        new_n = leaves[0];
        break;
      case 0x3u:
        new_n = ntk_dest.get_constant( true );
        new_n_const_or_ci = true;
        break;
      default:
        assert( false );
      }

      new_n_const_or_ci = ntk_dest.is_constant( ntk_dest.get_node( new_n ) ) || ntk_dest.is_ci( ntk_dest.get_node( new_n ) );

//      new_n_base_level = new_n_const_or_ci ? 0u : edge_level[{ ntk_dest.get_node( leaves[0] ), n }] + 1;
      new_n_base_level = new_n_const_or_ci ? 0u :  leaf_levels[0] + 1;
    }
    else
    {
      // should not have more than 4 fanin nodes
      assert( leaves.size() <= 4u );

      // if less than 4 fanins, add dummy inputs
      while ( leaves.size() < 4u )
      {
        leaves.push_back( ntk_dest.get_constant( false ) );
        leaf_levels.push_back( 0u );
        leaf_no_splitters.push_back( true );
      }

      auto tt = kitty::extend_to( f, 4u );

      // std::vector<uint32_t> fanin_levels( 4u, 0u );
      // std::vector<bool> fanin_const_or_ci( 4u, true );
      // for ( auto i = 0u; i < fanin.size(); i++ )
      // {
      //   fanin_const_or_ci[i] = ;
      //   fanin_levels[i] = ( fanin_const_or_ci[i] ? 0u : edge_level[{ ntk_dest.get_node( fanin[i] ), n }] );
      // }

      auto res = ( params.strategy == aqfp_node_resyn_strategy::best_cost )
                     ? db.get_best_cost_replacement( tt._bits[0], leaf_levels, leaf_no_splitters )
                     : db.get_best_level_replacement( tt._bits[0], leaf_levels, leaf_no_splitters );

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

        sigmap[x] = leaves[perm[ind]];
        levmap[x] = leaf_levels[perm[ind]];
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
        if ( ! (ntk_dest.is_constant( ntk_dest.get_node( sigmap[j] ) ) || ntk_dest.is_ci( ntk_dest.get_node( sigmap[j] )) )) {
        for ( auto k = 0u; k < repnet.nodes[j].size(); k++ )
        {
          auto lev_dif = ( gate_levels[repnet.nodes[j][k]] > gate_levels[j] ) ? gate_levels[repnet.nodes[j][k]] - gate_levels[j] : 1u;
          levmap[j] = std::max( levmap[j], levmap[repnet.nodes[j][k]] + lev_dif );
        }
        }

        // if we have already set the level, do not increase it
        level_update_callback(ntk_dest.get_node(sigmap[j]), levmap[j]);
      }

      new_n = output_inv ? !sigmap[0] : sigmap[0];
      new_n_const_or_ci = ntk_dest.is_constant( ntk_dest.get_node( new_n ) ) || ntk_dest.is_ci( ntk_dest.get_node( new_n ) );
      new_n_base_level = new_n_const_or_ci ? 0u : std::get<3>( res );
    }

    resyn_performed_callback( new_n);

  }

private:
  aqfp_node_resyn_param params;
  aqfp_logical_network_db<aqfp_logical_network_t<int>, 4u> db;
};

} // namespace mockturtle
