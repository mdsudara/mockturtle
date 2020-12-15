#pragma once

#include <chrono>
#include <fmt/format.h>
#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "./aqfp_cost_func.hpp"

using level_config_t = uint64_t;

void compute_all(
    std::vector<size_t>& level,
    const std::vector<size_t>& minlev,
    std::vector<size_t>& maxlev,
    const std::vector<std::vector<int>>& fanout,
    const mockturtle::aqfp_logical_network& net,
    size_t current_gid,
    double cost_so_far,
    std::unordered_map<level_config_t, double>& config_cost,
    size_t num_inputs = 4u )
{
  if ( net.nodes.size() == current_gid )
  {
    level_config_t lconfig = 0u;
    for ( auto f : net.input_slots )
    {
      if ( f == (int)net.zero_input )
      {
        continue;
      }

      lconfig = ( lconfig << 8 ) + level[f];
    }

    if ( !config_cost.count( lconfig ) )
    {
      config_cost[lconfig] = IMPOSSIBLE;
    }

    config_cost[lconfig] = std::min( config_cost[lconfig], cost_so_far );
    return;
  }

  if ( net.zero_input == current_gid )
  {
    compute_all( level, minlev, maxlev, fanout, net, current_gid + 1, cost_so_far, config_cost, num_inputs );
    return;
  }

  for ( auto lev = minlev[current_gid]; lev <= maxlev[current_gid]; lev++ )
  {
    auto cost = cost_for_node_if_in_level( lev, fanout[current_gid], level );
    if ( cost >= IMPOSSIBLE )
    {
      continue;
    }
    level[current_gid] = lev;
    compute_all( level, minlev, maxlev, fanout, net, current_gid + 1, cost_so_far + cost, config_cost, num_inputs );
  }
}

/**
 * Computes the AQFP cost for a network and save it if it is not computed before. Then returns
 * the computed cost. Cost is computed by considering different ways of fixing levels of the gates
 * and then computing locally optimum fanout-nets.
 */
void aqfp_cost_all( mockturtle::aqfp_logical_network& orig_net )
{
  assert( orig_net.is_dag );

  mockturtle::aqfp_logical_network net = orig_net;

  std::vector<size_t> minlev( net.nodes.size() );
  std::vector<size_t> maxlev( net.nodes.size() );
  std::vector<std::vector<int>> fanout( net.nodes.size() );

  // compute min levels from the top and fanouts
  compute_min_levels_and_fanouts( net, minlev, fanout );

  auto last_level = *( std::max_element( minlev.begin(), minlev.end() ) );

  std::vector<size_t> level( net.nodes.size(), 0u );

  std::unordered_map<level_config_t, double> config_cost;

  while ( true )
  {
    // compute max levels for gates from the bottom
    compute_max_levels( net, maxlev, fanout, last_level );

    // fix levels for the root and the inputs
    maxlev[0] = 0;

    // check different level configurations for the other gates and compute the cost for buffers and splitters
    compute_all( level, minlev, maxlev, fanout, net, 1u, 0.0, config_cost, 4 );

    if ( config_cost.size() > 0 )
    {
      break;
    }

    last_level++;
  }

  double cost_for_gates = compute_simple_cost( net, gate_costs );

  for ( auto it = config_cost.begin(); it != config_cost.end(); it++ )
  {
    fmt::print( "{:08x} {}\n", it->first, it->second + cost_for_gates );
  }
}
