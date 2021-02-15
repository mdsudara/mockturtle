#pragma once

#include <iostream>
#include <map>
#include <unordered_map>
#include <vector>

#include <kitty/kitty.hpp>

#include "./dag.hpp"
#include "./dag_cost.hpp"

/**
 *  TODO: [ ] For all real-primary-input configurations, keep entries in the database.
 *  TODO: [ ] Return levels of the gates together with the network so that we can build the real depth view.
 *  TODO: [ ] Return the correct inverter configuration so that we can construct the correct majority circuit. 
 */
namespace mockturtle
{
template<uint32_t N = 4u>
class npn_cache
{
  using npn_info = std::tuple<uint64_t, uint32_t, std::vector<uint8_t>>;

public:
  npn_cache() : arr( 1ul << ( 1ul << N ) ), has( 1ul << ( 1ul << N ), false )
  {
    static_assert( N == 4u, "Template parameter N must be 4 in the current implementation." );
  }

  npn_info operator()( uint64_t tt )
  {
    if ( has[tt] )
    {
      return arr[tt];
    }

    kitty::dynamic_truth_table dtt( N );
    dtt._bits[0] = tt;

    auto tmp = kitty::exact_npn_canonization( dtt );

    has[tt] = true;
    return ( arr[tt] = { std::get<0>( tmp )._bits[0] & 0xffff, std::get<1>( tmp ), std::get<2>( tmp ) } );
  }

private:
  std::vector<npn_info> arr;
  std::vector<bool> has;
};

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

template<typename Ntk, int N = 4>
class aqfp_logical_network_db
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

  aqfp_logical_network_db( double buffer_cost, uint32_t verbose = 0u )
      : buffer_cost( buffer_cost ), db( 1ul << ( 1ul << N ) ), verbose( verbose ), cc( { { 3u, 6.0 } }, { { 4u, 2.0 } }, buffer_cost, 4u )
  {
    static_assert( N == 4u, "Template parameter N must be 4 in the current implementation." );
  }

  /**
   * Filter database configurations that are "covered" by other configurations 
   */
  void remove_redundant()
  {
    for ( auto i = db.begin(); i != db.end(); i++ )
    {
      auto& npn = i->first;
      auto& configs = i->second;

      std::unordered_map<uint64_t, replacement> good_configs;
      for ( auto it = configs.begin(); it != configs.end(); it++ )
      {
        bool ok = true;
        for ( auto jt = configs.begin(); jt != configs.end(); jt++ )
        {
          if ( jt->first == it->first )
            continue;
          // check whether there exist input-wise smaller level config
          bool jt_smaller_to_it = true;
          for ( auto j = 0; j < N; j++ )
          {
            if ( level_of_input( jt->first, j ) > level_of_input( it->first, j ) )
            {
              jt_smaller_to_it = false;
              break;
            }
          }
          if ( jt_smaller_to_it )
          {
            double extra_buff_cost = 0;
            for ( auto j = 0; j < N; j++ )
            {
              extra_buff_cost += buffer_cost * ( level_of_input( it->first, j ) - level_of_input( jt->first, j ) );
            }
            if ( extra_buff_cost + jt->second.cost <= it->second.cost )
            {
              fmt::print( "[aqfp_db] configuration {:04x} already covered by {:04x} [{} {} {}].\n",
                          it->first, jt->first, it->second.cost, jt->second.cost, extra_buff_cost );
              ok = false;
              break;
            }
          }
        }
        if ( ok )
        {
          good_configs[it->first] = it->second;
        }
      }
      db[npn] = good_configs;
    }
  }

  /**
   * Update the database with a rusult for network `net`.
   */
  void update( const Ntk& ntk, const std::unordered_map<uint64_t, double>& cost_config )
  {
    std::vector<uint64_t> input_tt = { 0x0000UL, 0xaaaaUL, 0xccccUL, 0xf0f0UL, 0xff00UL };
    auto fs = all_functions_from_dag( input_tt, ntk );

    // for all functions synthesizable by assigning inputs in order, compute their npn class
    std::unordered_set<uint64_t> npn;
    for ( auto f : fs )
    {
      auto tmp = npndb( f );

      auto& npntt = std::get<0>( tmp );
      auto& npnperm = std::get<2>( tmp );

      if ( npn.count( npntt ) ) // already processed
      {
        continue;
      }
      npn.insert( npntt );

      // compute reverse mapping
      std::vector<uint8_t> revperm( N );
      for ( uint8_t i = 0; i < N; i++ )
      {
        revperm[npnperm[i]] = i;
      }

      for ( auto it = cost_config.begin(); it != cost_config.end(); it++ )
      {
        auto& lvl_cfg = it->first;
        auto& cost = it->second;
        // std::vector<uint8_t> new_levels = { level_of_input( lvl_cfg, revperm[0] ),
        //                                     level_of_input( lvl_cfg, revperm[1] ),
        //                                     level_of_input( lvl_cfg, revperm[2] ),
        //                                     level_of_input( lvl_cfg, revperm[3] ) };
        std::vector<uint8_t> new_levels = { level_of_input( lvl_cfg, npnperm[0] ),
                                            level_of_input( lvl_cfg, npnperm[1] ),
                                            level_of_input( lvl_cfg, npnperm[2] ),
                                            level_of_input( lvl_cfg, npnperm[3] ) };
        assert( new_levels.size() == N );
        auto new_lvl_cfg = lvl_cfg_from_vec( new_levels );

        if ( !db[npntt].count( new_lvl_cfg ) || db[npntt][new_lvl_cfg].cost > cost )
        {
          db[npntt][new_lvl_cfg] = { cost, ntk, new_levels, revperm, false, {}, {} };
        }
      }
    }
  }

  /**
   * Save database to the output stream `os`. 
   */
  void save_db_to_file( std::ostream& os )
  {
    // number of functions
    os << fmt::format( "{}\n", db.size() );
    for ( auto i = db.begin(); i != db.end(); i++ )
    {
      // npn class
      auto& k = i->first;

      auto& m = i->second;
      os << fmt::format( "{:04x}\n", k );

      // number of entries for the npn class
      os << fmt::format( "{}\n", m.size() );
      for ( auto it = m.begin(); it != m.end(); it++ )
      {
        auto& lvl_cfg = it->first;
        auto& r = it->second;

        //////////////////////
        // std::vector<uint8_t> perm( N );
        // for ( uint8_t i = 0; i < N; i++ )
        // {
        //   perm[r.input_perm[i]] = i;
        // }
        //  auto lvls = lvl_cfg_to_vec(lvl_cfg, 4u);
        //  std::vector<uint8_t> tmp(4u);
        //  for (auto j = 0u; j < 4u; j++) {
        //    tmp[j] = lvls[perm[j]];
        //  }
        //  std::reverse(tmp.begin(), tmp.begin() + r.ntk.input_slots.size() - (r.ntk.zero_input ? 1 : 0));
        //  for (auto j = 0u; j < 4u; j++) {
        //    lvls[j] = tmp[r.input_perm[j]];
        //  }
        //  auto lvl_cfg_new = lvl_cfg_from_vec(lvls);
        // os << fmt::format( "{:08x}\n", lvl_cfg_new);
        // UNCOMMENT NEXT LINE
        //////////////////////

        os << fmt::format( "{:08x}\n", lvl_cfg );
        os << fmt::format( "{}\n", r.cost );
        os << fmt::format( "{}\n", r.ntk.encode_as_string() );
        os << fmt::format( "{}\n", fmt::join( r.input_perm, " " ) );
      }
    }
  }

  /**
   * Load database from input stream `is`. 
   */
  void load_db_from_file( std::istream& is )
  {
    std::string line;

    std::getline( is, line );
    uint32_t num_func = std::stoul( line );

    for ( auto func = 0u; func < num_func; func++ )
    {
      std::getline( is, line );
      uint64_t npn = std::stoul( line, 0, 16 );

      std::getline( is, line );
      uint32_t num_entries = std::stoul( line );

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

        std::vector<uint8_t> perm( N );
        for ( int i = 0; i < N; i++ )
        {
          uint32_t t;
          is >> t;
          perm[i] = t;
        }
        std::getline( is, line ); // ignore the current line end

        // std::vector<uint8_t> perm( N );
        // for ( uint8_t i = 0; i < N; i++ )
        // {
        //   perm[r.input_perm[i]] = i;
        // }
        // auto lvls = lvl_cfg_to_vec(lvl_cfg, 4u);
        /////////////////////////////////
        // std::vector<uint8_t> tmp(4u);
        // for (auto j = 0u; j < 4u; j++) {
        //     tmp[j] = levels[perm[j]];
        // }
        // std::reverse(tmp.begin(), tmp.begin() + ntk.input_slots.size() - (ntk.zero_input ? 1 : 0));

        // for (auto j = 0u; j < 4u; j++) {
        //     levels[perm[j]] = tmp[j];
        // }
        /////////////////////////////////
        //  auto lvl_cfg_new = lvl_cfg_from_vec(lvls);
        // os << fmt::format( "{:08x}\n", lvl_cfg_new);
        // fix permutations
        // reverse permutation twice
        // std::vector<uint8_t> tmp(N);
        // for(int i = 0; i < N; i++) {
        //   tmp[perm[i]] = levels[i];
        // }

        //   std::reverse(tmp.begin(), tmp.end());
        //  for(int i = 0; i < N; i++) {
        //    levels[perm[i]] = tmp[i];
        //  }
        // levels = tmp;
        // std::cout << "reversing\n";
        ////////////////////////////

        lvl_cfg = lvl_cfg_from_vec( levels );

        if ( !db[npn].count( lvl_cfg ) || db[npn][lvl_cfg].cost > cost )
        {
          db[npn][lvl_cfg] = { cost, ntk, levels, perm, false, {}, {} };
        }
      }
    }
  }

  /**
   * returns (success, replacement, best cost, best_lev).
   */
  std::tuple<bool, replacement, double, uint32_t> get_best_cost_replacement(
      uint64_t f, std::vector<uint32_t> _levels, std::vector<bool> _is_const )
  {
    // find the npn class for the function
    auto tmp = npndb( f );
    auto& npntt = std::get<0>( tmp );
    auto& npnperm = std::get<2>( tmp );

    if ( db[npntt].empty() )
    {
      assert( false );
      return { false, replacement{}, 0.0, 0u };
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

      double cost = buffer_count * buffer_cost + r.cost;
      if ( cost < best_cost || ( cost == best_cost && max_lev < best_lev ) )
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

    auto [new_cost, new_levels] = cc.aqfp_cost_with_fixed_input_levels( best.ntk, levs, is_const );

    best.gate_levels = new_levels;
    best_cost += new_cost;

    return { true, fix_inverters_and_permutation( best, f ), best_cost, best_lev };
  }

  /**
   * returns (success, replacement, best cost, best_lev).
   */
  std::tuple<bool, replacement, double, uint32_t> get_best_level_replacement(
      uint64_t f, std::vector<uint32_t> _levels, std::vector<bool> _is_const )
  {
    // find the npn class for the function
    auto tmp = npndb( f );
    auto& npntt = std::get<0>( tmp );
    auto& npnperm = std::get<2>( tmp );

    if ( db[npntt].empty() )
    {
      assert( false );
      return { false, replacement{}, 0.0, 0u };
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
        if ( !is_const[i] )
        {
          auto temp = levels[i] + level_of_input( lvl_cfg, i );
          buffer_count += ( max_lev - temp );
        }
      }

      double cost = buffer_count * buffer_cost + r.cost;
      if ( max_lev < best_lev || ( max_lev == best_lev && cost < best_cost ) )
      {
        best_lev = max_lev;
        best_cost = cost;
        best = r;
      }
    }

    best_cost -= best.cost;

    std::vector<uint32_t> levs( 4u );
    for ( auto i = 0u; i < 4u; i++ )
    {
      levs[i] = best.input_levels[best.input_perm[i]];
    }

    auto [new_cost, new_levels] = cc.aqfp_cost_with_fixed_input_levels( best.ntk, levs, is_const );
    best.gate_levels = new_levels;
    best_cost += new_cost;

    return { true, fix_inverters_and_permutation( best, f ), best_cost, best_lev };
  }

  replacement fix_inverters_and_permutation( replacement rep, uint64_t func )
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

    return rep;
  }

private:
  double buffer_cost;
  std::unordered_map<uint64_t, std::unordered_map<uint64_t, replacement>> db;
  uint32_t verbose;
  npn_cache<N> npndb;
  mockturtle::aqfp_cost_computer<Ntk> cc;

  /** 
   * get all functions synthesizable from `net` if input slots are assigned the truthtables in `input_tt`.
   */
  std::unordered_set<uint64_t> all_functions_from_dag( const std::vector<uint64_t>& input_tt, const Ntk& net )
  {
    static_assert( N == 4u, "Template parameter N must be equal to 4 in the current implementation" );

    uint32_t num_inputs = net.input_slots.size();
    if ( net.zero_input != 0 )
    {
      num_inputs--;
    }

    std::unordered_set<uint64_t> res;

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
        tt[i - 1] = maj(
            ith_gate_config == 0ul ? ~tt[net.nodes[i - 1][0]] : tt[net.nodes[i - 1][0]],
            ith_gate_config == 1ul ? ~tt[net.nodes[i - 1][1]] : tt[net.nodes[i - 1][1]],
            ith_gate_config == 2ul ? ~tt[net.nodes[i - 1][2]] : tt[net.nodes[i - 1][2]] );
      }
      res.insert( tt[0] & 0xffff );
    }

    return res;
  }

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
        tt[i - 1] = maj(
            ith_gate_config == 0ul ? ~tt[net.nodes[i - 1][0]] : tt[net.nodes[i - 1][0]],
            ith_gate_config == 1ul ? ~tt[net.nodes[i - 1][1]] : tt[net.nodes[i - 1][1]],
            ith_gate_config == 2ul ? ~tt[net.nodes[i - 1][2]] : tt[net.nodes[i - 1][2]] );
      }

      // fmt::print("\tcomputed tt = {:04x}\n", tt[0] & 0xffff);

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

  /**
   * Compute the majority function on truth-tables.
   */
  template<typename T>
  inline static T maj( const T& a, const T& b, const T& c )
  {
    return ( a & b ) | ( c & ( a | b ) );
  }
};

} // namespace mockturtle
