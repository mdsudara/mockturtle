#pragma once

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include <kitty/kitty.hpp>

#include "./dag.hpp"
#include "./dag_cost.hpp"
#include "./gen_dag.hpp"
#include "./sat.hpp"
#include "./simulate_dag.hpp"

namespace mockturtle
{

template<uint32_t N = 4u>
class npn_cache
{
  using npn_info = std::tuple<uint64_t, uint32_t, std::vector<uint8_t>>;

public:
  npn_cache() : arr( 1ul << ( 1ul << N )), has( 1ul << ( 1ul << N ), false )
  {
  }

  uint64_t operator()( uint64_t tt )
  {
    if ( has[tt] )
      return arr[tt];
    kitty::dynamic_truth_table dtt( N );
    dtt._bits[0] = tt;
    has[tt] = true;
    return ( arr[tt] = kitty::exact_npn_canonization( dtt ) );
  }

private:
  std::vector<npn_info> arr;
  std::vector<bool> has;
};

/**
 * Returns the level of input with index `input_idx` from level configuration `lvl_cfg`.
 */
uint8_t level_of_input( uint64_t lvl_cfg, uint32_t input_idx )
{
  return ( lvl_cfg >> ( 8u * input_idx ) ) & 0xff;
}

template<typename Ntk>
struct replacement
{
  double cost;
  Ntk ntk;
  std::vector<uint8_t> inverters;
  std::vector<uint8_t> levels;
};

template<typename Ntk, int N = 4>
class aqfp_logical_network_db
{
  using lvl_cfg_t = uint64_t;

public:
  aqfp_logical_network_db( uint32_t verbose = 0u ) : verbose( verbose ), db(1ul << (1ul << N))
  {
    assert( N == 4u );

    std::vector<uint8_t> indices( N );
    std::iota( indices.begin(), indices.end(), 0u );

    std::vector<uint64_t> input_tt = { 0x0000UL, 0xaaaaUL, 0xccccUL, 0xf0f0UL, 0xff00UL };
    do
    {
      perm.push_back( indices );
      input_perm.push_back( input_tt );
    } while ( std::next_permutation( indices.begin(), indices.end() ) && std::next_permutation( input_tt.begin() + 1, input_tt.end() ) );
    assert( perm.size() == 24u );
    assert( input_perm.size() == 24u );
  }

  void update( const Ntk& ntk, const std::unordered_map<uint64_t, double>& cost_config )
  {
    for ( auto i = 0u; i < perm.size(); i++ )
    {
      auto fs = all_functions_from_dag( input_perm[i], ntk );
      for ( auto f : fs )
      {
        f >>= 1;
        for ( auto it = cost_config.begin(); it != cost_config.end(); it++ )
        {
          auto& lvl_cfg = it->first;
          auto& cost = it->second;
          std::vector<uint8_t> new_levels = { level_of_input( lvl_cfg, perm[i][0] ),
                                              level_of_input( lvl_cfg, perm[i][1] ),
                                              level_of_input( lvl_cfg, perm[i][2] ),
                                              level_of_input( lvl_cfg, perm[i][3] ) };
          assert( new_levels.size() == N );
          auto new_lvl_cfg = lvl_cfg_from_vec( new_levels );
          if ( !db[f].count( new_lvl_cfg ) || db[f][new_lvl_cfg].cost > cost )
          {
            Ntk ntk2 = ntk;

            // begin permute inputs
            decltype( ntk.input_slots ) input_slots;
            ntk2.input_slots.clear();
            for ( auto j = 0u; j < ntk.input_slots.size(); j++ )
            {
              if ( ntk.input_slots[j] != ntk.zero_input )
              {
                input_slots.push_back( ntk.input_slots[j] );
              }
            }
            for ( auto j = 0u; j < N; j++ )
            {
              ntk2.input_slots.push_back( input_slots[perm[i][j]] );
            }
            if ( ntk.zero_input != 0 )
            {
              ntk2.input_slots.push_back( ntk.zero_input );
            }
            // end permute inputs

            db[f][new_lvl_cfg] = { cost, ntk2, {}, new_levels };
          }
        }
      }
    }
  }

  /** 
   * Get all functions synthesizable from `net` if input slots are assigned the truthtables in `input_tt`.
   */
  std::unordered_set<uint64_t> all_functions_from_dag( const std::vector<uint64_t>& input_tt, const Ntk& net )
  {
    uint32_t num_inputs = net.input_slots.size();
    if ( net.zero_input != 0 )
    {
      num_inputs--;
    }

    std::unordered_set<uint64_t> res;

    std::vector<uint64_t> tt( net.nodes.size(), input_tt[0] );
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

  /**
   * Save the databse to a file.
   */
  void save_db_to_file( std::ostream& os )
  {
    // number of functions
    os << fmt::format( "{}\n", db.size() );
    for ( auto& m : db )
    {
      // number of entries per function
      os << fmt::format( "{}\n", m.size() );
      for ( auto it = m.begin(); it != m.end(); it++ )
      {
        auto& lvl_cfg = it->first;
        auto& r = it->second;
        os << fmt::format( "{:08x}\n", lvl_cfg );
        os << fmt::format( "{}\n", r.cost );
        os << fmt::format( "{}\n", r.ntk.encode_as_string() );
      }
    }
  }

  /**
   * Load the database from a file.
   */
  void load_db_from_file( std::istream& is )
  {
    uint32_t num_func;
    is >> num_func;
    for ( auto func = 0u; func < num_func; func++ )
    {
      uint64_t num_entries;
      is >> num_entries;

      std::string line;
      for ( auto j = 0u; j < num_entries; j++ )
      {
        std::getline( is, line );
        uint64_t lvl_cfg = std::stoul( line, 0, 16 );
        auto levels = lvl_cfg_to_vec( lvl_cfg );

        std::getline( is, line );
        double cost = std::stod( line );

        std::getline( is, line );
        Ntk ntk;
        ntk.decode_dag( line );
        db[func][lvl_cfg] = { cost, ntk, {}, levels };
      }
    }
  }

  /**
   * Returns (success, replacement, buffer cost).
   */
  std::tuple<bool, replacement<Ntk>, double> get_best_replacement( uint64_t f, std::vector<uint32_t> levels, double buffer_cost )
  {
    assert( (f & 1) == 0u );
    f >>= 1;

    if ( db[f].empty() )
      return { false, replacement<Ntk>{}, 0.0 };

    double best_cost = std::numeric_limits<double>::infinity();
    replacement<Ntk>& best = db[f].begin()->second;

    for ( auto it = db[f].begin(); it != db[f].end(); it++ )
    {
      auto& lvl_cfg = it->first;
      auto& r = it->second;

      uint32_t max_lev = 0u;
      for ( auto i = 0u; i < N; i++ )
      {
        auto temp = levels[i] + level_of_input( lvl_cfg, i );
        max_lev = std::max( max_lev, temp );
      }
      uint32_t buffer_count = 0u;
      for ( auto i = 0u; i < N; i++ )
      {
        auto temp = levels[i] + level_of_input( lvl_cfg, i );
        buffer_count += ( max_lev - temp );
      }

      double cost = buffer_count * buffer_cost + r.cost;
      if ( cost < best_cost )
      {
        best_cost = cost;
        best = r;
      }
    }

    return { true, best, best_cost };
  }

private:
  uint32_t verbose;
  std::vector<std::map<lvl_cfg_t, replacement<Ntk>>> db;
  std::vector<std::vector<uint8_t>> perm;
  std::vector<std::vector<uint64_t>> input_perm;

  static std::vector<uint8_t> lvl_cfg_to_vec( uint64_t lvl_cfg )
  {
    std::vector<uint8_t> res( N );
    for ( auto i = 0u; i < N; i++ )
    {
      res[i] = level_of_input( lvl_cfg, i );
    }
    return res;
  }

  static uint64_t lvl_cfg_from_vec( std::vector<uint8_t> levels )
  {
    uint64_t res = 0u;
    for ( auto i = 0u; i < levels.size(); i++ )
    {
      res = ( res << 8u ) + levels[i];
    }
    return res;
  }
};

} // namespace mockturtle
