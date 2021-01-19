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
  npn_cache() : arr( 1ul << ( 1ul << N ) ), has( 1ul << ( 1ul << N ), false )
  {
    assert( N == 4u );
  }

  npn_info operator()( uint64_t tt )
  {
    if ( has[tt] )
    {
      return arr[tt];
    }

    _size++;

    kitty::dynamic_truth_table dtt( N );
    dtt._bits[0] = tt;

    auto tmp = kitty::exact_npn_canonization( dtt );

    has[tt] = true;
    return ( arr[tt] = { std::get<0>( tmp )._bits[0] & 0xffff, std::get<1>( tmp ), std::get<2>( tmp ) } );
  }

  uint32_t size()
  {
    return _size;
  }

private:
  std::vector<npn_info> arr;
  std::vector<bool> has;
  uint32_t _size = 0u;
};

/**
 * Returns the level of input with index `input_idx` from level configuration `lvl_cfg`.
 */
uint8_t level_of_input( uint64_t lvl_cfg, uint32_t input_idx )
{
  return ( lvl_cfg >> ( 8u * input_idx ) ) & 0xff;
}

/**
 * Returns the vector representation of the level configuration `lvl_cfg`.
 */
std::vector<uint8_t> lvl_cfg_to_vec( uint64_t lvl_cfg, uint32_t num_leaves )
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
uint64_t lvl_cfg_from_vec( std::vector<uint8_t> levels )
{
  uint64_t res = 0u;
  for ( auto i = 0u; i < levels.size(); i++ )
  {
    res |= ( levels[i] << ( 8u * i ) );
  }
  return res;
}

template<typename Ntk>
struct replacement
{
  double cost;
  Ntk ntk;
  std::vector<uint8_t> inverters;
  std::vector<uint8_t> levels;
  std::vector<uint8_t> perm;
};

template<typename Ntk, int N = 4>
class aqfp_logical_network_db
{
public:
  aqfp_logical_network_db( uint32_t verbose = 0u ) : verbose( verbose ), db( 1ul << ( 1ul << N ) )
  {
    assert( N == 4u );
  }

  void update( const Ntk& ntk, const std::unordered_map<uint64_t, double>& cost_config )
  {
    std::vector<uint64_t> input_tt = { 0x0000UL, 0xaaaaUL, 0xccccUL, 0xf0f0UL, 0xff00UL };
    auto fs = all_functions_from_dag( input_tt, ntk );

    std::unordered_set<uint64_t> npn;

    // for all functions synthesizable by assigning inputs in order, compute their npn class
    for ( auto f : fs )
    {
      auto tmp = npndb( f );

      auto& npntt = std::get<0>( tmp );
      auto& npnperm = std::get<2>( tmp );

      if ( npn.count( npntt ) )
      {
        continue;
      }

      npn.insert( npntt );

      // compute reverse mapping
      std::vector<uint8_t> q( N );
      for ( uint8_t i = 0; i < N; i++ )
      {
        q[npnperm[i]] = i;
      }

      for ( auto it = cost_config.begin(); it != cost_config.end(); it++ )
      {
        auto& lvl_cfg = it->first;
        auto& cost = it->second;
        std::vector<uint8_t> new_levels = { level_of_input( lvl_cfg, q[0] ),
                                            level_of_input( lvl_cfg, q[1] ),
                                            level_of_input( lvl_cfg, q[2] ),
                                            level_of_input( lvl_cfg, q[3] ) };
        assert( new_levels.size() == N );
        auto new_lvl_cfg = lvl_cfg_from_vec( new_levels );

        if ( !db[npntt].count( new_lvl_cfg ) || db[npntt][new_lvl_cfg].cost > cost )
        {
          db[npntt][new_lvl_cfg] = { cost, ntk, {}, new_levels, q };
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

    auto xx = net.input_slots;
    std::sort( xx.begin(), xx.end() );
    assert( xx == net.input_slots );

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
      res.insert( tt[0] & 0xffff );
    }

    return res;
  }

  /**
   * Save the database to a file.
   */
  void save_db_to_file( std::ostream& os )
  {
    // number of functions
    os << fmt::format( "{}\n", db.size() );
    for ( auto i = db.begin(); i != db.end(); i++ )
    {
      // npn class
      auto& k = i->first;
      os << fmt::format( "{:04x}\n", k );

      auto& m = i->second;
      // number of entries for the npn class
      os << fmt::format( "{}\n", m.size() );
      for ( auto it = m.begin(); it != m.end(); it++ )
      {
        auto& lvl_cfg = it->first;
        auto& r = it->second;
        os << fmt::format( "{:08x}\n", lvl_cfg );
        os << fmt::format( "{}\n", r.cost );
        os << fmt::format( "{}\n", r.ntk.encode_as_string() );
        os << fmt::format( "{}\n", fmt::join( r.perm, " " ) );
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
      uint64_t npn;
      is >> npn;

      uint32_t num_entries;
      is >> num_entries;

      std::string line;
      std::getline( is, line ); // ignore the current line end
      for ( auto j = 0u; j < num_entries; j++ )
      {
        std::getline( is, line );
        uint64_t lvl_cfg = std::stoul( line, 0, 16 );
        auto levels = lvl_cfg_to_vec( lvl_cfg, N );

        std::getline( is, line );
        double cost = std::stod( line );

        std::getline( is, line );
        Ntk ntk;
        ntk.decode_dag( line );

        std::vector<uint8_t> perm;
        for ( int i = 0; i < N; i++ )
        {
          is >> perm[i];
        }

        if ( !db[npn].count( lvl_cfg ) || db[npn][lvl_cfg].cost > cost )
        {
          db[npn][lvl_cfg] = { cost, ntk, {}, levels, perm };
        }
        std::cerr << "db npn class count " << db.size() << "\n";
      }
    }
  }

  /** 
   * Return the network after permuting the input slots of `net` according to `perm`.
   */
  Ntk get_input_permuted_ntk( const Ntk& net, const std::vector<uint8_t>& perm )
  {
    Ntk res = net;
    res.input_slots.clear();

    decltype( net.input_slots ) input_slots;
    for ( auto j = 0u; j < net.input_slots.size(); j++ )
    {
      if ( net.input_slots[j] != net.zero_input )
      {
        input_slots.push_back( net.input_slots[j] );
      }
    }
    for ( auto j = 0u; j < N; j++ )
    {
      res.input_slots.push_back( ( perm[j] < input_slots.size() ) ? input_slots[perm[j]] : net.zero_input );
    }
    if ( net.zero_input != 0 )
    {
      res.input_slots.push_back( net.zero_input );
    }

    std::cerr << fmt::format( "{} | {}\n", fmt::join( input_slots, " " ), fmt::join( res.input_slots, " " ) );
    std::cerr << fmt::format( "{} | {}\n", fmt::join( input_slots, " " ), fmt::join( net.input_slots, " " ) );
    std::cerr << "dbg " << res.encode_as_string() << "\n";
    return res;
  }

  /**
   * Returns (success, replacement, buffer cost).
   */
  std::tuple<bool, replacement<Ntk>, double> get_best_replacement( uint64_t f, std::vector<uint32_t> _levels, double buffer_cost )
  {
    auto tmp = npndb( f );
    auto& npntt = std::get<0>( tmp );
    auto& npnperm = std::get<2>( tmp );

    if ( db[npntt].empty() )
      return { false, replacement<Ntk>{}, 0.0 };

    std::vector<uint32_t> levels( _levels.size() );

    for ( int i = 0; i < levels.size(); i++ )
    {
      levels[i] = _levels[npnperm[i]];
    }

    double best_cost = std::numeric_limits<double>::infinity();
    replacement<Ntk>& best = db[f].begin()->second;

    for ( auto it = db[npndb].begin(); it != db[npndb].end(); it++ )
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
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, replacement<Ntk>>> db;
  npn_cache<N> npndb;

}; // namespace mockturtle

} // namespace mockturtle
