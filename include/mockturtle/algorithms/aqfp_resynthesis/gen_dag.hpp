#pragma once

#include <chrono>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <vector>

#include <fmt/format.h>

#include "./dag.hpp"
#include "./dag_builder.hpp"
#include "./gen_dag_util.hpp"

namespace mockturtle
{

using part = std::multiset<int>;
using partition = std::multiset<part>;
using partition_set = std::set<partition>;

struct dag_generator_params
{

  uint32_t max_gates;      // max number of gates allowed
  uint32_t max_num_fanout; // max number of fanouts per gate
  uint32_t max_width;      // max number of gates in any given level
  uint32_t max_num_in;     // max number of primary input slots (including the constant)
  uint32_t max_levels;     // max number of gate levels

  std::vector<uint32_t> allowed_num_fanins;        // the types of allowed majority gates
  std::map<uint32_t, uint32_t> max_gates_of_fanin; // max number of gates allowed for each type

  bool allow_par_costing; // allow parallelized cost computations
  uint32_t cost_threads;
  uint32_t verbose;

  dag_generator_params() : max_gates( std::numeric_limits<uint32_t>::max() ),
                           max_num_fanout( std::numeric_limits<uint32_t>::max() ),
                           max_width( std::numeric_limits<uint32_t>::max() ),
                           max_num_in( std::numeric_limits<uint32_t>::max() ),
                           allowed_num_fanins( { 3u } ),
                           max_gates_of_fanin( { { 3u, std::numeric_limits<uint32_t>::max() } } ),
                           allow_par_costing( false ),
                           cost_threads( 1u ),
                           verbose( 0u ) {}
};

template<typename CostFn>
class dag_compare
{
public:
  dag_compare( CostFn cc ) : cc( cc )
  {
  }

  /*! \brief Returns true if s should appear before f.
   *
   * This happens if s has a smaller cost or if s has the same cost as f
   * but s is a DAG whereas f is partial DAG.
   */
  template<typename T>
  bool operator()( T& f, T& s )
  {
    return ( f.second > s.second ) || ( f.second == s.second && !s.first.is_partial_dag && f.first.is_partial_dag );
  }

private:
  CostFn cc;
};

template<typename NodeT, typename CostFn>
class dag_generator
{
  using NtkBuilder = aqfp_dag_builder<NodeT>;
  using Ntk = aqfp_dag<NodeT>;

public:
  dag_generator( const dag_generator_params& params, CostFn cc ) : params( params ), cc( cc ), pq( dag_compare( cc ) )
  {
    for ( auto&& fin : params.allowed_num_fanins )
    {
      if ( params.max_gates_of_fanin.at( fin ) > 0 )
      {
        auto root = NtkBuilder::get_root( fin );
        pq.push( { root, cc( root ) } );
      }
    }
  }

  template<typename PredicateT>
  std::optional<std::pair<Ntk, double>> next_dag( PredicateT&& should_expand )
  {
    while ( true )
    {
      if ( pq.empty() )
      {
        return std::nullopt;
      }

      auto [res, res_cost] = pq.top();
      pq.pop();

      if ( !res.is_partial_dag )
      {
        return {{ res, res_cost} };
      }

      if ( params.max_levels > res.num_levels )
      {
        auto ext = get_layer_extension( res );
        std::vector<double> costs( ext.size() );

        if ( params.allow_par_costing )
        {
          auto block_size = params.cost_threads;

          for ( auto i = 0u; i < ext.size(); i += block_size )
          {
            std::vector<std::future<double>> temp;
            for ( auto j = i; j < i + block_size && j < ext.size(); j++ )
            {
              temp.push_back( std::async(
                  std::launch::async, [&]( NtkBuilder& item ) { return cc( item ); }, std::ref( ext[j] ) ) );
            }

            for ( auto j = i; j < i + block_size && j < ext.size(); j++ )
            {
              costs[j] = temp[j - i].get();
            }
          }
        }

        for ( auto i = 0u; i < ext.size(); i++ )
        {
          pq.push( { ext[i], costs[i] } );
        }
      }

      if ( should_expand( res ) )
      {
        auto dags = get_dags_from_partial_dag( res );
        std::vector<double> costs(dags.size());

        if ( params.allow_par_costing )
        {
          auto block_size = params.cost_threads;

          for ( auto i = 0u; i < dags.size(); i += block_size )
          {
            std::vector<std::future<double>> temp;
            for ( auto j = i; j < i + block_size && j < dags.size(); j++ )
            {
              temp.push_back( std::async( std::launch::async, [&]( NtkBuilder& item ) { return cc( item ); }, std::ref( dags[j] ) ) );
            }

            for ( auto j = i; j < i + block_size && j < dags.size(); j++ )
            {
              costs[j] = temp[j - i].get();
            }
          }
        }

        for ( auto i = 0u; i < dags.size(); i++ )
        {
          pq.push( {dags[i], costs[i]} );
        }
      }

      if ( params.verbose > 0u )
      {
        std::cerr << fmt::format( "expanded partial dag {}\ncurrent size of pq {}\n", as_string( res ), pq.size() );
      }
    }
  }

  std::optional<std::pair<Ntk, double>> next_dag()
  {
    static auto always_expand = []( auto& net ) { (void) net; return true; };
    return next_dag( always_expand );
  }

  /*! \brief Extend the current aqfp logical network by one more level. */
  std::vector<NtkBuilder> get_layer_extension( const NtkBuilder& net )
  {
    // Choose which last layer slots to use
    // Choose which non-last later slots to use
    // get partitions of the chosen last layer slots
    // extend the partitions using chosen non-last layer slots

    std::vector<NtkBuilder> result;

    auto max_counts = net.max_equal_fanins();

    auto last_options = sublist_gen( net.last_layer_leaves );
    auto other_options = sublist_gen( net.other_leaves );

    const auto last_counts = detail::get_frequencies( net.last_layer_leaves );
    const auto other_counts = detail::get_frequencies( net.other_leaves );

    /* Consider all different ways of choosing a non-empty subset of last layer slots */
    for ( auto&& last : last_options )
    {
      if ( last.empty() )
        continue;

      /* Remaining slots in the last layer */
      auto last_counts_cpy = last_counts;
      for ( auto&& e : last )
      {
        last_counts_cpy[e]--;
      }

      /* Consider all different ways of choosing a subset of other layer slots */
      for ( auto&& other : other_options )
      {

        /* Remaining slots in the other layers */
        auto other_counts_cpy = other_counts;
        for ( auto&& e : other )
        {
          other_counts_cpy[e]--;
        }

        /* Compute the new set of other leaves for all resulting partial DAGs */
        std::vector<int> other_leaves_new;
        for ( auto it = last_counts_cpy.begin(); it != last_counts_cpy.end(); it++ )
        {
          for ( auto i = 0u; i < it->second; i++ )
          {
            other_leaves_new.push_back( it->first );
          }
        }
        for ( auto it = other_counts_cpy.begin(); it != other_counts_cpy.end(); it++ )
        {
          for ( auto i = 0u; i < it->second; i++ )
          {
            other_leaves_new.push_back( it->first );
          }
        }

        if ( params.max_gates == 0 || params.max_gates > net.num_gates() )
        {
          auto max_gates = params.max_gates > 0u ? params.max_gates - net.num_gates() : 0u;
          auto last_layers_partitions = partition_gen( last, max_counts, max_gates, params.max_num_fanout );

          for ( auto p : last_layers_partitions )
          {
            auto extensions = partition_ext( other, p, max_counts, params.max_num_fanout );
            for ( auto q : extensions )
            {
              auto temp = get_next_partial_dags( net, q, other_leaves_new );
              for ( auto&& r : temp )
              {
                r.num_levels++;
                result.push_back( r );
              }
            }
          }
        }
      }
    }

    return result;
  }

  /**
   * @brief Compute all DAGs derived from a given partial DAG.
   */
  std::vector<NtkBuilder> get_dags_from_partial_dag( const NtkBuilder& net )
  {
    std::vector<NodeT> leaves = net.last_layer_leaves;
    leaves.insert( leaves.end(), net.other_leaves.begin(), net.other_leaves.end() );

    std::sort( leaves.begin(), leaves.end() );

    auto max_counts = net.max_equal_fanins();

    auto partitions = partition_gen( leaves, max_counts, params.max_num_in, 0 /* unlimited part sizes (fanouts) */ );

    std::vector<NtkBuilder> result;
    for ( auto p : partitions )
    {
      auto new_net = get_next_dag( net, p );
      result.push_back( new_net );

      for ( auto i = 0u; i < new_net.input_slots.size(); i++ )
      {
        auto temp_net = new_net;
        temp_net.zero_input = new_net.input_slots[i];
        result.push_back( temp_net );
      }
    }

    return result;
  }

  /**
   * \brief Compute the partial DAGs obtained by combining the slots of 'orig' as
   * indicated by partitioning 'p'. We have to try different gate types for each part
   * in partition 'p'.
   */
  std::vector<NtkBuilder> get_next_partial_dags( const NtkBuilder& orig, const partition& p, const std::vector<int>& other_leaves )
  {
    auto max_allowed_of_fanin = params.max_gates_of_fanin;
    for ( auto it = orig.num_gates_of_fanin.begin(); it != orig.num_gates_of_fanin.end(); it++ )
    {
      assert( max_allowed_of_fanin[it->first] >= it->second );
      max_allowed_of_fanin[it->first] -= it->second;
    }

    std::vector<part> q( p.begin(), p.end() );

    /* For each part in partition 'p', consider different types of gates to connect. */
    auto res = add_node_recur( orig, q, 0, max_allowed_of_fanin );

    for ( auto&& net : res )
    {
      net.other_leaves = other_leaves;
    }

    return res;
  }

  /**
   * \brief Recursively consider different fanin gates for different parts in partition 'p'.
   */
  std::vector<NtkBuilder> add_node_recur( const NtkBuilder& orig, const std::vector<part>& p, uint32_t ind, std::map<uint32_t, uint32_t>& max_allowed_of_fanin )
  {
    if ( ind == p.size() )
    {
      return { orig.copy_without_leaves() };
    }

    std::vector<NtkBuilder> res;

    /* Decide what fanin gate to use for part in partition 'p' at index 'ind'. */
    for ( auto&& fin : params.allowed_num_fanins )
    {
      if ( max_allowed_of_fanin[fin] == 0 )
      {
        continue;
      }
      max_allowed_of_fanin[fin]--;

      auto temp = add_node_recur( orig, p, ind + 1, max_allowed_of_fanin );

      for ( auto&& t : temp )
      {
        auto net = t.copy_with_last_layer_leaves();
        net.add_internal_node( fin, p[ind], true );
        res.push_back( net );
      }

      max_allowed_of_fanin[fin]++;
    }

    return res;
  }

  /**
   * \brief Compute the DAG obtained from partial DAG 'orig' by combining slots
   * according to partitions 'p'.
   */
  NtkBuilder get_next_dag( const NtkBuilder& orig, const partition& p )
  {
    assert( orig.is_partial_dag );

    auto net = orig.copy_without_leaves();
    net.is_partial_dag = false;
    std::for_each( p.begin(), p.end(), [&net]( auto q ) { net.add_leaf_node( q ); } );

    return net;
  }

private:
  uint32_t next_pdag_id = 0u;

  dag_generator_params params;
  CostFn cc;
  std::priority_queue<std::pair<NtkBuilder, double>, std::vector<std::pair<NtkBuilder, double>>, dag_compare<CostFn>> pq;

  detail::partition_generator<NodeT> partition_gen;
  detail::partition_extender<NodeT> partition_ext;
  detail::sublist_generator<NodeT> sublist_gen;
};

} // namespace mockturtle
