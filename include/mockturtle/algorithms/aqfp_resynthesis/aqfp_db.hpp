#pragma once

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include <kitty/kitty.hpp>

#include "./dag.hpp"
#include "./dag_cost.hpp"
#include "./npn_cache.hpp"

namespace mockturtle
{

template<typename T>
T bitwise_majority( const T& a, const T& b, const T& c )
{
  return ( a & b ) | ( c & ( a | b ) );
}

template<typename T>
T bitwise_majority( const T& a, const T& b, const T& c, const T& d, const T& e )
{
  return ( a & b & c ) | ( a & b & d ) | ( a & b & e ) | ( a & c & d ) | ( a & c & e ) |
         ( a & d & e ) | ( b & c & d ) | ( b & c & e ) | ( b & d & e ) | ( c & d & e );
}

/**
 * Returns the level of input with index `input_idx` from level configuration `lvl_cfg`.
 */
inline uint8_t level_of_input( uint64_t lvl_cfg, uint32_t input_idx )
{
  return ( lvl_cfg >> ( 8u * input_idx ) ) & 0xff;
}

/**
 * Returns the vector representation of the level configuration `lvl_cfg`.
 */
inline std::vector<uint8_t> lvl_cfg_to_vec( uint64_t lvl_cfg, uint32_t num_leaves )
{
  std::vector<uint8_t> res( num_leaves );
  for ( auto i = 0u; i < num_leaves; i++ )
  {
    res[i] = level_of_input( lvl_cfg, i );
  }
  return res;
}

/**
 * Returns the level configuration for levels represented by `levels`.
 */
inline uint64_t lvl_cfg_from_vec( std::vector<uint8_t> levels )
{
  uint64_t res = 0u;
  for ( auto i = 0u; i < levels.size(); i++ )
  {
    res |= ( levels[i] << ( 8u * i ) );
  }
  return res;
}

template<typename Ntk = aqfp_dag<>, int N = 4>
class aqfp_db
{

public:
  struct replacement
  {
    double cost;
    Ntk ntk;
    std::vector<uint8_t> input_levels;    // input levels
    std::vector<uint8_t> input_perm;      // input permutation so that ntk compute the npn-class
    bool output_inv;                      // output invertedness so that the ntk computes the npn-class
    std::vector<uint8_t> inverter_config; // inverter config of the gates so that the ntk computs the npn-class
    std::vector<uint32_t> gate_levels;    // levels of the internal gates
  };

  aqfp_db(
      const std::unordered_map<uint32_t, double>& gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } },
      const std::unordered_map<uint32_t, double>& splitters = { { 1u, 2.0 }, { 4u, 2.0 } } )
      : gate_costs( gate_costs ), splitters( splitters ), cc( gate_costs, splitters, 4u )
  {
    static_assert( N == 4u, "Template parameter N must be 4 in the current implementation." );
  }

  aqfp_db( const std::unordered_map<uint64_t, std::map<uint64_t, replacement>>& db,
           const std::unordered_map<uint32_t, double>& gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } },
           const std::unordered_map<uint32_t, double>& splitters = { { 1u, 2.0 }, { 4u, 2.0 } } )
      : gate_costs( gate_costs ), splitters( splitters ), db( db ), cc( gate_costs, splitters )
  {
    static_assert( N == 4u, "Template parameter N must be 4 in the current implementation." );
  }

  using gate_info = std::vector<uint32_t>;                                               // fanin list with lsb denoting the inversion
  using mig_structure = std::tuple<std::vector<gate_info>, std::vector<uint32_t>, bool>; // (gates, levels, output inverted flag);

  /**
   * returns (success, replacement, best cost, best_lev).
   */
  template<typename ComparisonFn>
  mig_structure get_best_replacement( uint64_t f, std::vector<uint32_t> _levels,
                                      std::vector<bool> _is_const, ComparisonFn&& comparison_fn )
  {
    // find the npn class for the function
    auto tmp = npndb( f );
    auto& npntt = std::get<0>( tmp );
    auto& npnperm = std::get<2>( tmp );

    if ( db[npntt].empty() )
    {
      assert( false );
      return { {}, {}, false };
    }

    // map input levels
    std::vector<uint32_t> levels( _levels.size() );
    std::vector<bool> is_const( _levels.size() );
    for ( auto i = 0u; i < levels.size(); i++ )
    {
      levels[i] = _levels[npnperm[i]];
      is_const[i] = _is_const[npnperm[i]];
    }

    double best_cost = std::numeric_limits<double>::infinity();
    uint32_t best_lev = std::numeric_limits<uint32_t>::max();
    replacement best = db[npntt].begin()->second;

    for ( auto it = db[npntt].begin(); it != db[npntt].end(); it++ )
    {
      const auto& lvl_cfg = it->first;
      const auto& r = it->second;

      uint32_t max_lev = 0u;
      for ( auto i = 0u; i < N; i++ )
      {
        auto temp = levels[i] + level_of_input( lvl_cfg, i );
        max_lev = std::max( max_lev, temp );
      }
      uint32_t buffer_count = 0u;
      for ( auto i = 0u; i < N; i++ )
      {
        if ( !is_const[i] )
        {
          auto temp = levels[i] + level_of_input( lvl_cfg, i );
          buffer_count += ( max_lev - temp );
        }
      }

      double cost = buffer_count * splitters.at( 1u ) + r.cost;

      if ( comparison_fn( { cost, max_lev }, { best_cost, best_lev } ) )
      {
        best_cost = cost;
        best_lev = max_lev;
        best = r;
      }
    }

    best_cost -= best.cost;

    std::vector<uint32_t> levs( 4u );
    for ( auto i = 0u; i < 4u; i++ )
    {
      levs[i] = best.input_levels[best.input_perm[i]];
    }

    auto [new_cost, new_levels] = cc( best.ntk, levs  );

    best.gate_levels = new_levels;
    best_cost += new_cost;

    return fix_inverters_and_permutation( best, f );
  }

private:
  std::unordered_map<uint32_t, double> gate_costs;
  std::unordered_map<uint32_t, double> splitters;
  std::unordered_map<uint64_t, std::map<uint64_t, replacement>> db;
  mockturtle::dag_aqfp_cost_and_depths<Ntk> cc;
  npn_cache<N> npndb;

  std::vector<uint8_t> inverter_config_for_func( const std::vector<uint64_t>& input_tt, const Ntk& net, uint64_t func )
  {
    uint32_t num_inputs = net.input_slots.size();
    if ( net.zero_input != 0 )
    {
      num_inputs--;
    }

    std::vector<uint64_t> tt( net.nodes.size(), input_tt[0] );
    auto input_ind = 1u;

    auto tmp_input_slots = net.input_slots;
    std::sort( tmp_input_slots.begin(), tmp_input_slots.end() );
    assert( tmp_input_slots == net.input_slots );

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
        tt[i - 1] = bitwise_majority(
            ith_gate_config == 0ul ? ~tt[net.nodes[i - 1][0]] : tt[net.nodes[i - 1][0]],
            ith_gate_config == 1ul ? ~tt[net.nodes[i - 1][1]] : tt[net.nodes[i - 1][1]],
            ith_gate_config == 2ul ? ~tt[net.nodes[i - 1][2]] : tt[net.nodes[i - 1][2]] );
      }

      if ( func == ( tt[0] & 0xffff ) )
      {
        std::vector<uint8_t> res;
        for ( auto i = 0u; i < net.num_gates(); i++ )
        {
          res.push_back( ( inv_config >> ( 2 * i ) ) & 3u );
        }
        return res;
      }
    }

    assert( false );
    return {};
  }

  mig_structure fix_inverters_and_permutation( replacement rep, uint64_t func )
  {
    auto [npntt, npn_inv, npnperm] = npndb( func );

    std::vector<uint8_t> ind = { 0u, 1u, 2u, 3u };

    std::vector<uint8_t> ind_func_from_npn = {
        ind[npnperm[0]],
        ind[npnperm[1]],
        ind[npnperm[2]],
        ind[npnperm[3]] };

    std::vector<uint8_t> ind_func_from_dag = {
        ind_func_from_npn[rep.input_perm[0]],
        ind_func_from_npn[rep.input_perm[1]],
        ind_func_from_npn[rep.input_perm[2]],
        ind_func_from_npn[rep.input_perm[3]] };

    rep.input_perm = ind_func_from_dag;

    std::vector<uint64_t> input_tt = {
        0xaaaaUL,
        0xccccUL,
        0xf0f0UL,
        0xff00UL };

    std::vector<uint64_t> input_perm_tt = {
        0x0000UL,
        input_tt[rep.input_perm[0]],
        input_tt[rep.input_perm[1]],
        input_tt[rep.input_perm[2]],
        input_tt[rep.input_perm[3]] };

    if ( func & 1 )
    {
      rep.output_inv = true;
      func = ( ~func ) & 0xffff;
    }

    rep.inverter_config = inverter_config_for_func( input_perm_tt, rep.ntk, func );

    std::vector<gate_info> gates{ {}, {}, {}, {}, {} };
    std::vector<uint32_t> depths{ 0u, 0u, 0u, 0u, 0u };

    std::map<uint32_t, uint32_t> sigmap;

    std::vector<uint32_t> inputs;
    auto i = 0u;
    for ( auto x : rep.ntk.input_slots )
    {
      if ( x == rep.ntk.zero_input )
      {
        sigmap[x] = 0u;
      }
      else
      {
        sigmap[x] = rep.input_perm[i++] + 1;
      }
      depths[sigmap[x]] = rep.gate_levels[x];
    }

    for ( auto i = rep.ntk.num_gates(); i > 0; i-- )
    {
      auto j = i - 1;
      sigmap[j] = rep.ntk.num_gates() + 4u - j;

      gates.push_back( {} );
      assert( gates.size() == sigmap[j] + 1 );
      depths.push_back( rep.gate_levels[j] );

      auto type = rep.inverter_config[j];
      auto& node = rep.ntk.nodes[j];
      for ( auto k = 0u; k < node.size(); k++ )
      {
        auto new_fanin_id = sigmap[node[k]];
        auto new_fanin_inv = is_kth_fanin_inverted( node.size(), k, type );
        gates[sigmap[j]].push_back( ( new_fanin_id << 1 ) | ( new_fanin_inv ? 1u : 0u ) );
      }
    }

    return { gates, depths, rep.output_inv };
  }

  static bool is_kth_fanin_inverted( uint32_t num_fanin, uint32_t fanin_idx, uint32_t type )
  {
    if ( num_fanin == 3u )
    {
      return type == fanin_idx;
    }
    else
    {
      assert( false );
      return false;
    }
  }
};

} // namespace mockturtle
