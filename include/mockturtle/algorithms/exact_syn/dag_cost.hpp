#pragma once

#include <numeric>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>

#include "./dag.hpp"
#include "./gen_dag_util.hpp"

namespace mockturtle
{

template<typename Ntk>
class simple_cost_computer
{
public:
  simple_cost_computer( const std::unordered_map<uint32_t, double>& gate_costs ) : gate_costs( gate_costs )
  {
  }

  /**
   * \brief Computes and updates the simple cost of a network.
   */
  double operator()( mockturtle::aqfp_logical_network_t<int>& net )
  {
	if ( net.computed_cost < 0.0 )
	{
	  net.computed_cost = compute_simple_cost( net );
	}
	return net.computed_cost;
  }

private:
  const std::unordered_map<uint32_t, double> gate_costs;

  double compute_simple_cost( const Ntk& net )
  {
	double res = 0.0;

	for ( auto it = net.num_gates_of_fanin.begin(); it != net.num_gates_of_fanin.end(); it++ )
	{
	  assert( gate_costs.count( it->first ) > 0 );
	  res += gate_costs.at( it->first ) * it->second;
	}

	return res;
  }
};

template<typename Ntk>
class aqfp_cost_computer
{
  using level_config_t = uint64_t;

  static constexpr double IMPOSSIBLE = std::numeric_limits<double>::infinity();

public:
  aqfp_cost_computer(
	  const std::unordered_map<uint32_t, double>& gate_costs,
	  const std::vector<std::pair<uint32_t, double>>& splitters,
	  double buffer_cost,
	  uint32_t max_num_pis ) : simp_cc( gate_costs ), splitters( splitters ), buffer_cost( buffer_cost ), max_num_pis( max_num_pis ) {}

  double operator()( Ntk& net )
  {
	return aqfp_cost( net );
  }

  /**
   * Computes the AQFP cost for a network and save it if it is not computed before. Then returns
   * the computed cost. Cost is computed by considering different ways of fixing levels of the gates
   * and then computing locally optimum fanout-nets.
   */
  double aqfp_cost( Ntk& orig_net )
  {
	if ( orig_net.computed_cost >= 0.0 )
	{
	  return orig_net.computed_cost;
	}

	cost_context ctx( orig_net.is_partial_dag ? convert_to_min_cost_dag( orig_net ) : orig_net );

	compute_fanouts( ctx.net, ctx.fanout );
	compute_min_levels( ctx.net, ctx.fanout, ctx.minlev );

	auto lastlev = *( std::max_element( ctx.minlev.begin(), ctx.minlev.end() ) );

	double cost = IMPOSSIBLE;
	while ( true )
	{
	  compute_max_levels( ctx.net, ctx.fanout, ctx.maxlev, lastlev );

	  // fix levels for the root and the inputs
	  ctx.maxlev[0] = 0;
	  for ( auto f : ctx.net.input_slots )
	  {
		if ( f != ctx.net.zero_input )
		{
		  ctx.minlev[f] = lastlev;
		}
	  }

	  // check different level configurations for the other gates and compute the cost for buffers and splitters
	  cost = compute_best_cost( ctx, 1u, 0.0 );

	  if ( cost < IMPOSSIBLE )
	  {
		break;
	  }

	  lastlev++;
	}

	// add the gate costs
	cost += simp_cc( ctx.net );

	return ( orig_net.computed_cost = cost );
  }

  /**
   * Computes the AQFP cost for a network and save it if it is not computed before. Then returns
   * the computed cost. Cost is computed by considering different ways of fixing levels of the gates
   * and then computing locally optimum fanout-nets.
   */
  std::unordered_map<level_config_t, double>
  cost_all_level_configurations( const Ntk& orig_net )
  {
	assert( !orig_net.is_partial_dag );

	std::unordered_map<level_config_t, double> config_cost;

	cost_context ctx( orig_net );

	compute_fanouts( ctx.net, ctx.fanout );
	compute_min_levels( ctx.net, ctx.fanout, ctx.minlev );

	auto lastlev = *( std::max_element( ctx.minlev.begin(), ctx.minlev.end() ) );
	while ( true )
	{
	  compute_max_levels( ctx.net, ctx.fanout, ctx.maxlev, lastlev );

	  // fix levels for the root and the inputs
	  ctx.maxlev[0] = 0;

	  // check different level configurations for the other gates and compute the cost for buffers and splitters
	  compute_best_costs_for_all_configs( ctx, 1u, 0.0, config_cost );

	  if ( config_cost.size() > 0 )
	  {
		break;
	  }

	  lastlev++;
	}

	double cost_for_gates = simp_cc( ctx.net );

	for ( auto it = config_cost.begin(); it != config_cost.end(); it++ )
	{
	  it->second += cost_for_gates;
	}

	return config_cost;
  }

private:
  struct cost_context
  {
	Ntk net;
	std::vector<uint32_t> curlev;
	std::vector<uint32_t> minlev;
	std::vector<uint32_t> maxlev;
	std::vector<std::vector<typename Ntk::NodeT>> fanout;
	std::unordered_map<std::vector<uint32_t>, double> rellev_cache;

	cost_context( const Ntk& net ) : net( net ), curlev( net.nodes.size() ), minlev( net.nodes.size() ), maxlev( net.nodes.size() ), fanout( net.nodes.size() ) {}
  };

  simple_cost_computer<Ntk> simp_cc;
  const std::vector<std::pair<uint32_t, double>> splitters;
  double buffer_cost;
  uint32_t max_num_pis;

  /**
   * \brief Computes the min-cost dag that can be obtained from the given partial dag.
   */
  Ntk convert_to_min_cost_dag( Ntk net )
  {
	assert( net.is_partial_dag );

	// we need to compute the cost of a min cost dag that can be obtained from this partial dag
	// for this, we ignore up to max number of fanins that could be connected to constants, and
	// connect remaining slots to distinct inputs

	auto mef = net.max_equal_fanins();

	std::multiset<typename Ntk::NodeT> zero = {};
	for ( auto&& t : net.last_layer_leaves )
	{
	  if ( mef[t] > 0 )
	  {
		mef[t]--;
		zero.insert( t );
	  }
	  else
	  {
		net.add_leaf_node( { t } );
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
		net.add_leaf_node( { t } );
	  }
	}

	net.zero_input = net.add_leaf_node( zero );

	net.last_layer_leaves.clear();
	net.other_leaves.clear();

	return net;
  }

  /**
   * \brief Computes the fanouts of nodes.
   */
  template<typename FanOutT>
  void compute_fanouts( const Ntk& net, FanOutT& fanout )
  {
	for ( auto i = 0u; i < net.nodes.size(); i++ )
	{
	  for ( auto&& f : net.nodes[i] )
	  {
		if ( net.zero_input != f )
		{
		  fanout[f].push_back( i );
		}
	  }
	}
  }

  /**
   * \brief Computes the min levels of nodes.
   */
  template<typename FanOutT, typename MinLevT>
  void compute_min_levels( const Ntk& net, const FanOutT& fanout, MinLevT& minlev )
  {
	for ( auto i = 0u; i < net.nodes.size(); i++ )
	{
	  if ( fanout[i].size() == 0u )
	  {
		minlev[i] = 0u;
	  }
	  else
	  {
		auto critical_fo = *( std::max_element( fanout[i].begin(), fanout[i].end(),
												[&minlev]( auto x, auto y ) { return ( minlev[x] < minlev[y] ); } ) );
		minlev[i] = 1 + minlev[critical_fo];
		if ( fanout[i].size() > 1 )
		{
		  minlev[i]++;
		}
	  }
	}
  }

  /**
   * \brief Computes the max levels from bottom assuming the input slots are at 'lastlev'.
   */
  template<typename FanOutT, typename MaxLevT>
  void compute_max_levels( const Ntk& net, const FanOutT& fanout, MaxLevT& maxlev, uint32_t lastlev )
  {
	for ( auto& f : net.input_slots )
	{
	  if ( f != net.zero_input )
	  {
		maxlev[f] = lastlev;
	  }
	  else
	  {
		maxlev[f] = std::numeric_limits<uint32_t>::max();
	  }
	}

	for ( auto i = net.num_gates(); i > 0u; i-- )
	{
	  maxlev[i - 1] = std::numeric_limits<uint32_t>::max();
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

  /**
   * \brief Compute the best splitter and buffer cost for a given relative level configuration 'config'.
   */
  double cost_for_config( cost_context& ctx, const std::vector<uint32_t> config )
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

	if ( ctx.rellev_cache.count( config ) )
	{
	  return ctx.rellev_cache[config];
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

		temp += cost_for_config( ctx, new_config );

		if ( temp < result )
		{
		  result = temp;
		}
	  }
	}

	return ( ctx.rellev_cache[config] = result );
  }

  /**
   * \brief Find the locally optimal fanout-net for a node if it in level 'lev', if its fanouts are
   * 'fanouts', and if fanouts are placed at fixed level given by 'level'.
   */
  double cost_for_node_if_in_level( cost_context& ctx, uint32_t lev, std::vector<typename Ntk::NodeT> fanouts )
  {
	std::vector<uint32_t> rellev;
	for ( auto fo : fanouts )
	{
	  rellev.push_back( lev - ctx.curlev[fo] );
	}
	std::sort( rellev.begin(), rellev.end() );

	return cost_for_config( ctx, rellev );
  }

  /**
   * \brief Compute the AQFP cost for a given network. This function fixes the level of gate with index 'current_gid',
   * and calls it recursively for the next gate id. Levels are fixed starting from the root in the BFS order,
   * so once we fix a level of a gate, we can compute the cost for its fanout-net. Once all gate levels are
   * fixed, we can compute the costs of fanout-nets for the inputs (unless the network is partial DAG).
   * For partial DAGs, we still add dummy-inputs to make them look like DAGs, but we ignore the costs of their
   * fanout-nets.
   */
  double compute_best_cost( cost_context& ctx, uint32_t current_gid, double cost_so_far )
  {
	if ( ctx.net.num_gates() == current_gid )
	{
	  if ( ctx.net.is_partial_dag )
	  {
		auto t = ctx.net.input_slots.size() - max_num_pis + 1;

		// let k be the least number of splitters we need
		// k (sp.size() - 1) + 1 >= t ===> k >= (t - 1) / (sp.size() - 1)

		// we cannot say anything about the buffer requirement not knowning the actual inputs

		return cost_so_far + ( ( t - 1 ) / ( splitters.back().first - 1 ) ) * splitters.back().second;
	  }

	  // net is a DAG, so compute fanout-net costs for primary inputs

	  for ( auto f : ctx.net.input_slots )
	  {
		if ( f == ctx.net.zero_input )
		{
		  continue;
		}

		cost_so_far += cost_for_node_if_in_level( ctx, ctx.maxlev[f], ctx.fanout[f] );

		if ( cost_so_far >= IMPOSSIBLE )
		{
		  return IMPOSSIBLE;
		}
	  }

	  return cost_so_far;
	}

	auto result = IMPOSSIBLE;
	for ( auto lev = ctx.minlev[current_gid]; lev <= ctx.maxlev[current_gid]; lev++ )
	{
	  auto cost = cost_for_node_if_in_level( ctx, lev, ctx.fanout[current_gid] );
	  if ( cost >= IMPOSSIBLE )
	  {
		continue;
	  }
	  ctx.curlev[current_gid] = lev;
	  auto temp = compute_best_cost( ctx, current_gid + 1, cost_so_far + cost );
	  if ( temp < result )
	  {
		result = temp;
	  }
	}

	return result;
  }

  /**
   * \brief Compute costs considering all input level configurations. Works for at most 8 inputs and at most 255 levels.
   */
  void compute_best_costs_for_all_configs( cost_context& ctx, uint32_t current_gid, double cost_so_far, std::unordered_map<level_config_t, double>& config_cost )
  {
	if ( ctx.net.nodes.size() == current_gid )
	{
	  level_config_t lconfig = 0u;
	  for ( auto f : ctx.net.input_slots )
	  {
		if ( f == ctx.net.zero_input )
		{
		  continue;
		}

		lconfig = ( lconfig << 8 ) + ctx.curlev[f];
	  }

	  if ( !config_cost.count( lconfig ) )
	  {
		config_cost[lconfig] = IMPOSSIBLE;
	  }

	  config_cost[lconfig] = std::min( config_cost[lconfig], cost_so_far );
	  return;
	}

	if ( ctx.net.zero_input == (typename Ntk::NodeT)current_gid )
	{
	  compute_best_costs_for_all_configs( ctx, current_gid + 1, cost_so_far, config_cost );
	  return;
	}

	for ( auto lev = ctx.minlev[current_gid]; lev <= ctx.maxlev[current_gid]; lev++ )
	{
	  auto cost = cost_for_node_if_in_level( ctx, lev, ctx.fanout[current_gid] );
	  if ( cost >= IMPOSSIBLE )
	  {
		continue;
	  }
	  ctx.curlev[current_gid] = lev;
	  compute_best_costs_for_all_configs( ctx, current_gid + 1, cost_so_far + cost, config_cost );
	}
  }
};

} // namespace mockturtle

