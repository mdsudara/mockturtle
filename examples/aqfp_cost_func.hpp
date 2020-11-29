#pragma once

#include <chrono>
#include <fmt/format.h>
#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <numeric>
#include <vector>

const double IMPOSSIBLE = 1.0e+20;

// costs for gates of different number of fanins
std::map<size_t, double> gate_costs = {{3u, 3.0}, {5u, 5.0}};

// avialable splitter fanouts and thier costs
std::vector<std::pair<size_t, double>> splitters = {{3u, 3.0}};

// cost of a buffer
double buffer_cost = 1.0;

/**
 * Compute the cost of network considering only the gate costs.
 */
double compute_simple_cost( const mockturtle::aqfp_logical_network& net, const std::map<size_t, double>& gates )
{
  double res = 0.0;

  for ( auto it = net.num_gates_of_fanin.begin(); it != net.num_gates_of_fanin.end(); it++ )
  {
    if ( gates.count( it->first ) )
    {
      res += gates.at( it->first ) * it->second;
    }
  }

  return res;
}

/**
 * Computes and updates the simple cost of a network.
 */
double simple_cost( mockturtle::aqfp_logical_network& net )
{
  if ( net.computed_cost < 0.0 )
  {
    net.computed_cost = compute_simple_cost( net, gate_costs );
  }

  return net.computed_cost;
}

/**
 * Returns the maximum value of an item in 'items'.
 */
template<typename T, typename V>
size_t max_val_of_items( T&& items, const V& values )
{
  size_t res = 0u;

  for ( auto&& t : items )
  {
    if ( values[t] > res )
    {
      res = values[t];
    }
  }

  return res;
}

/**
 * Returns the minimum value of an item in 'items'.
 */
template<typename T, typename V>
size_t min_val_of_items( T&& items, const V& values )
{
  size_t res = 1000000000u;

  for ( auto&& t : items )
  {
    if ( values[t] < res )
    {
      res = values[t];
    }
  }

  return res;
}

/**
 * Find and return the max element of a part in parition 'p'.
 */
int max_in_partition( const mockturtle::partition& p )
{
  int res = -1;

  for ( const auto& q : p )
  {
    for ( const auto& r : q )
    {
      if ( r > res )
      {
        res = r;
      }
    }
  }

  return res;
}

/**
 * Computes the min-cost dag that can be obtained from the given partial dag.
 */
void convert_to_min_cost_dag( mockturtle::aqfp_logical_network& net )
{
  assert( net.is_partial_dag );

  // we need to compute the cost of a min cost dag that can be obtained from this partial dag
  // for this, we ignore up to max number of fanins that could be connected to constants, and
  // connect remaining slots to distinct inputs

  auto mef = net.max_equal_fanins();

  mockturtle::part zero = {};
  for ( auto&& t : net.last_layer_leaves )
  {
    if ( mef[t] > 0 )
    {
      mef[t]--;
      zero.insert( t );
    }
    else
    {
      net.add_leaf_node( {t} );
    }
  }

  for ( auto&& t : net.other_leaves )
  {
    if ( mef[t] > 0 )
    {
      mef[t]--;
      zero.insert( t );
    }
    else
    {
      net.add_leaf_node( {t} );
    }
  }

  net.zero_input = net.add_leaf_node( zero );

  net.last_layer_leaves.clear();
  net.other_leaves.clear();
}

/**
 * Computes the min levels (from top), fanouts, and number of real fanins (ignoring constants) of nodes in a network.
 */
void compute_min_levels_and_fanouts( const mockturtle::aqfp_logical_network& net, std::vector<size_t>& minlev, std::vector<std::vector<int>>& fanout )
{
  for ( auto i = 0u; i < net.nodes.size(); i++ )
  {
    for ( auto&& f : net.nodes[i] )
    {
      if ( (int)net.zero_input != f )
      {
        fanout[f].push_back( i );
      }
    }

    if ( fanout[i].size() == 0u )
    {
      minlev[i] = 0u;
    }
    else
    {
      minlev[i] = 1 + max_val_of_items( fanout[i], minlev );
      if ( fanout[i].size() > 1 )
      {
        minlev[i]++;
      }
    }
  }
}

void compute_max_levels( const mockturtle::aqfp_logical_network& net, std::vector<size_t>& maxlev, std::vector<std::vector<int>>& fanout, size_t last_lev )
{
  for ( auto& f : net.input_slots )
  {
    if ( f != (int)net.zero_input )
    {
      maxlev[f] = last_lev;
    }
    else
    {
      maxlev[f] = 1000000000u;
    }
  }

  for ( auto i = net.num_gates(); i > 0u; i-- )
  {
    maxlev[i - 1] = 1000000000u;
    for ( auto f : net.nodes[i - 1] )
    {
      auto t = maxlev[f] - 1;
      if ( fanout[f].size() > 1 )
      {
        t--;
      }
      if ( t < maxlev[i - 1] )
      {
        maxlev[i - 1] = t;
      }
    }
  }
}

thread_local std::map<std::vector<size_t>, double> config_cache;

double cost_for_config( const std::vector<size_t> config )
{
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

  if ( !config_cache.count( config ) )
  {
    auto result = IMPOSSIBLE;

    for ( const auto& s : splitters )
    {
      for ( auto size = 2u; size <= std::min( s.first, config.size() ); size++ )
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

        std::vector<size_t> new_config( config.begin(), config.begin() + ( config.size() - size ) );
        new_config.push_back( sp_lev );
        std::sort( new_config.begin(), new_config.end() );

        temp += cost_for_config( new_config );

        if ( temp < result )
        {
          result = temp;
        }
      }
    }

    config_cache[config] = result;
  }

  return config_cache[config];
}

/**
 * Find the locally optimal fanout-net for a node if it in level 'lev', if its fanouts are 
 * 'fanouts', and if fanouts are placed at fixed level given by 'level'.
 */
template<typename FanOutT, typename LevelT>
double cost_for_node_if_in_level( size_t lev, const FanOutT& fanouts, LevelT& level )
{
  std::vector<size_t> rellev;
  for ( auto fo : fanouts )
  {
    rellev.push_back( lev - level[fo] );
  }

  std::sort( rellev.begin(), rellev.end() );

  return cost_for_config( rellev );
}

/**
 * Compute the AQFP cost for a given network. This function fixes the level of gate with index 'current_gid', 
 * and calls it recursively for the next gate id. Levels are fixed starting from the root in the BFS order, 
 * so once we fix a level of a gate, we can compute the cost for its fanout-net. Once all gate levels are 
 * fixed, we can compute the costs of fanout-nets for the inputs (unless the network is partial DAG).
 * For partial DAGs, we still add dummy-inputs to make them look like DAGs, but we ignore the costs of their
 * fanout-nets.
 */
double compute(
    std::vector<size_t>& level,
    const std::vector<size_t>& minlev,
    std::vector<size_t>& maxlev,
    const std::vector<std::vector<int>>& fanout,
    const mockturtle::aqfp_logical_network& net,
    size_t current_gid,
    size_t last_gid,
    double cost_so_far,
    size_t num_inputs = 4u )
{
  if ( last_gid == current_gid )
  {
    if ( net.is_partial_dag )
    {
      auto t = net.input_slots.size() - num_inputs + 1;

      // let k be the least number of splitters we need
      // k (sp.size() - 1) + 1 >= t ===> k >= (t - 1) / (sp.size() - 1)

      // we cannot say anything about the buffer requirement not knowning the actual inputs

      return cost_so_far + ( ( t - 1 ) / ( splitters.back().first - 1 ) ) * splitters.back().second;
    }
    else
    {
      for ( auto f : net.input_slots )
      {
        if ( f == (int)net.zero_input )
        {
          continue;
        }

        cost_so_far += cost_for_node_if_in_level( maxlev[f], fanout[f], level );

        if ( cost_so_far >= IMPOSSIBLE )
        {
          return IMPOSSIBLE;
        }
      }
    }

    return cost_so_far;
  }

  auto result = IMPOSSIBLE;
  for ( auto lev = minlev[current_gid]; lev <= maxlev[current_gid]; lev++ )
  {
    auto cost = cost_for_node_if_in_level( lev, fanout[current_gid], level );
    if ( cost >= IMPOSSIBLE )
    {
      continue;
    }
    level[current_gid] = lev;
    auto temp = compute( level, minlev, maxlev, fanout, net, current_gid + 1, last_gid, cost_so_far + cost, num_inputs );
    if ( temp < result )
    {
      result = temp;
    }
  }

  return result;
}

/**
 * Computes the AQFP cost for a network and save it if it is not computed before. Then returns
 * the computed cost. Cost is computed by considering different ways of fixing levels of the gates
 * and then computing locally optimum fanout-nets.
 */
double aqfp_cost( mockturtle::aqfp_logical_network& orig_net )
{
  if ( orig_net.computed_cost < 0.0 )
  {
    mockturtle::aqfp_logical_network net = orig_net;
    if ( net.is_partial_dag )
    {
      convert_to_min_cost_dag( net );
    }

    std::vector<size_t> minlev( net.nodes.size() );
    std::vector<size_t> maxlev( net.nodes.size() );
    std::vector<std::vector<int>> fanout( net.nodes.size() );

    // compute min levels from the top and fanouts
    compute_min_levels_and_fanouts( net, minlev, fanout );

    auto last_level = *( std::max_element( minlev.begin(), minlev.end() ) );

    std::vector<size_t> level( net.nodes.size(), 0u );

    double cost = IMPOSSIBLE;
    while ( true )
    {
      // compute max levels for gates from the bottom
      compute_max_levels( net, maxlev, fanout, last_level );

      // fix levels for the root and the inputs
      maxlev[0] = 0;
      for ( auto f : net.input_slots )
      {
        if ( f != (int)net.zero_input )
        {
          minlev[f] = last_level;
        }
      }

      // check different level configurations for the other gates and compute the cost for buffers and splitters
      cost = compute( level, minlev, maxlev, fanout, net, 1u, net.num_gates(), 0.0 );

      if ( cost < IMPOSSIBLE )
      {
        break;
      }

      last_level++;
    }

    // add the gate costs
    cost += compute_simple_cost( net, gate_costs );

    orig_net.computed_cost = cost;
  }

  return orig_net.computed_cost;
}
