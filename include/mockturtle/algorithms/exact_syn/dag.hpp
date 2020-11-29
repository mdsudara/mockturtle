#pragma once

#include <chrono>
#include <fmt/format.h>
#include <future>
#include <map>
#include <mutex>
#include <queue>
#include <set>
#include <vector>

namespace mockturtle
{

using part = std::multiset<int>;
using partition = std::multiset<part>;
using partition_set = std::set<partition>;

/*
 * Computes and returns the string representation of a partition of elements
 */
std::string as_string( const partition& p )
{
  std::vector<std::string> parts;
  for ( auto&& q : p )
  {
    parts.push_back( fmt::format( "[{}]", fmt::join( q, ", " ) ) );
  }
  return fmt::format( "[ {} ]", fmt::join( parts, " " ) );
}

namespace detail
{

/*
 * Computes and returns the frequency map for a given collection of elements
 */
auto get_frequencies( const std::vector<int>& elems )
{
  std::map<int, size_t> elem_counts;
  for ( const auto& e : elems )
  {
    elem_counts[e]++;
  }
  return elem_counts;
}

using partition_cache_key_t = std::tuple<const std::vector<int>, const std::map<int, size_t>, size_t, size_t>;
thread_local std::map<partition_cache_key_t, partition_set> partition_cache;
/*
 * Computes and returns a vector of partitions for a given list of elements
 * such that no part contains any element 'e' more than 'max_counts[e]' times. 
 * Ex: if 'elems' = [ 1, 2, 2, 3 ], this will generate the partitions,
 * [ [1] [2] [2] [3] ], [ [1] [2] [2, 3] ], [ [1, 2] [2] [3] ],
 * [ [1, 2], [2, 3] ], [ [1, 2, 3], [2] ], and [ [1, 3] [2] [2] ], if
 * 'max_counts' is all-ones.
 */
partition_set get_all_partitions(
    std::vector<int> elems,
    const std::map<int, size_t>& max_counts = {},
    size_t max_parts = 0,
    size_t max_part_size = 0 )
{
  if ( elems.size() == 0 )
  {
    return { {} }; // return the empty parition.
  }

  partition_cache_key_t key = { elems, max_counts, max_parts, max_part_size };
  if ( !partition_cache.count( key ) )
  {
    partition_set result;

    auto last = elems.back();
    elems.pop_back();

    auto temp = get_all_partitions( elems, max_counts, max_parts, max_part_size );

    for ( auto&& t : temp )
    {
      partition cpy;

      // take 'last' in its own partition
      cpy = t;

      if ( max_parts == 0 || max_parts > cpy.size() )
      {
        cpy.insert( { last } );
        result.insert( cpy );
      }

      // add 'last' to one of the existing partitions
      for ( auto it = t.begin(); it != t.end(); )
      {
        if ( !max_counts.count( last ) || it->count( last ) < max_counts.at( last ) )
        {

          if ( max_part_size == 0 || max_part_size > it->size() )
          {
            cpy = t;
            auto elem_it = cpy.find( *it );
            auto cpy_elem = *elem_it;
            cpy_elem.insert( last );
            cpy.erase( elem_it );
            cpy.insert( cpy_elem );
            result.insert( cpy );
          }
        }

        std::advance( it, t.count( *it ) );
      }
    }

    partition_cache[key] = result;
  }

  return partition_cache[key];
}

using partition_ext_cache_key_t = std::tuple<const std::vector<int>, partition, const std::map<int, size_t>, size_t>;
thread_local std::map<partition_ext_cache_key_t, partition_set> partition_ext_cache;
/*
 * Compute a list of different partitions that can be obtained by adding elements in
 * 'elems' to the parts of 'base' such that no part contains any element 'e' more than
 * 'max_counts[e]' times.
 * Note: Assumes 'base' does not share elements with 'elems'.
 * Ex: if 'base' = [ [1, 2] [3] ] and 'elems' = [4, 4, 5], this will generate the
 * partitoins, [ [1, 2, 4] [3, 4, 5] ] and [ [1, 2, 4, 5] [3, 4] ], if
 * 'max_counts' is all-ones.
 */
partition_set extend_partitions( std::vector<int> elems, partition base, const std::map<int, size_t>& max_counts, size_t max_part_size = 0 )
{
  if ( elems.size() == 0 )
  {
    return { base };
  }

  partition_ext_cache_key_t key = { elems, base, max_counts, max_part_size };
  if ( !partition_ext_cache.count( key ) )
  {
    partition_set result;

    auto last = elems.back();
    elems.pop_back();

    auto temp = extend_partitions( elems, base, max_counts, max_part_size );
    for ( auto&& t : temp )
    {
      partition cpy;

      for ( auto it = t.begin(); it != t.end(); )
      {
        if ( it->count( last ) < max_counts.at( last ) )
        {

          if ( max_part_size == 0 || max_part_size > it->size() )
          {
            cpy = t;
            auto elem_it = cpy.find( *it );
            auto cpy_elem = *elem_it;
            cpy_elem.insert( last );
            cpy.erase( elem_it );
            cpy.insert( cpy_elem );
            result.insert( cpy );
          }
        }

        std::advance( it, t.count( *it ) );
      }
    }

    partition_ext_cache[key];
  }

  return partition_ext_cache[key];
}

using sub_list_cache_key_t = std::map<int, size_t>;
thread_local std::map<sub_list_cache_key_t, std::set<std::vector<int>>> sub_list_cache;
/*
 * Given element frequencies 'elem_counts', generate all "sub_lists" containing each element 'e' at
 * most 'elem_counts[e]' times.
 */
std::set<std::vector<int>> get_sub_lists_recur( std::map<int, size_t> elem_counts )
{
  if ( elem_counts.size() == 0u )
  {
    return { {} };
  }

  sub_list_cache_key_t key = elem_counts;
  if ( !sub_list_cache.count( key ) )
  {
    auto last = std::prev( elem_counts.end() );
    auto last_elem = last->first;
    auto last_count = last->second;
    elem_counts.erase( last );

    std::set<std::vector<int>> result;

    std::vector<int> t;
    for ( auto i = last_count; i > 0; --i )
    {
      t.push_back( last_elem );
      result.insert( t ); // insert a copy of t, and note that t is already sorted.
    }

    auto temp = get_sub_lists_recur( elem_counts );

    for ( std::vector<int> t : temp )
    {
      result.insert( t );
      for ( auto i = last_count; i > 0; --i )
      {
        t.push_back( last_elem );
        std::sort( t.begin(), t.end() );
        result.insert( t );
      }
    }

    sub_list_cache[key] = result;
  }

  return sub_list_cache[key];
}

/*
 * Given a list of elements 'elems', generate all sub lists of those elements.
 * Ex: if 'elems' = [1, 2, 2, 3], this will generate the following lists:
 * [0], [1], [1, 2], [1, 2, 2], [1, 2, 2, 3], [1, 2, 3], [1, 3], [2], [2, 2], [2, 2, 3], [2, 3], and [3].
 */
std::set<std::vector<int>> get_sub_lists( std::vector<int> elems )
{
  auto elem_counts = get_frequencies( elems );
  return get_sub_lists_recur( elem_counts );
}

} // namespace detail

struct dag_params
{

  size_t max_gates;
  size_t max_num_fanout;
  size_t max_width;
  size_t max_num_in;
  size_t max_level;
  std::vector<size_t> allowed_num_fanins;
  std::map<size_t, size_t> max_gates_of_fanin;

  dag_params() : max_gates( 1000u ), max_num_fanout( 1000u ), max_width( 1000u ), max_num_in( 1000u ), allowed_num_fanins( { 3u, 5u } ), max_gates_of_fanin( { { 3u, 1000u }, { 5u, 1000u } } )
  {
  }

  dag_params( size_t max_gates, size_t max_num_fanout, size_t max_width, size_t max_num_in )
      : max_gates( max_gates ), max_num_fanout( max_num_fanout ), max_width( max_width ), max_num_in( max_num_in )
  {
  }
};

/*
 * Enocdes the logical structure of an aqfp network.
 */
template<typename node_t = int>
struct aqfp_logical_network_t
{
  size_t pdag_id = 0u;

  bool is_partial_dag = false;
  bool is_dag = false;

  double computed_cost = -1.0;
  size_t level = 0u;

  std::vector<std::vector<node_t>> nodes;
  std::vector<size_t> node_num_fanin;
  std::map<size_t, size_t> num_gates_of_fanin;

  std::vector<node_t> last_layer_leaves;
  std::vector<node_t> other_leaves;
  std::vector<node_t> input_slots;

  size_t zero_input = 0;

  /*
   * Return the number of majority gates.
   */
  size_t num_gates() const
  {
    return nodes.size() - input_slots.size();
  }

  std::map<node_t, size_t> max_equal_fanins() const
  {
    std::map<node_t, size_t> res;
    for ( auto i = 0u; i < num_gates(); i++ )
    {
      res[i] = node_num_fanin[i] / 2;
    }
    return res;
  }

  void add_fanin( node_t node_id, node_t fanin )
  {
    nodes[node_id].push_back( fanin );
  }

  // add a new node with given number of fanins
  size_t add_internal_node( size_t num_fanin = 3u, const part& fanouts = {}, bool is_in_last_layer = true )
  {
    size_t node_id = nodes.size();
    nodes.push_back( {} );

    assert( node_id == node_num_fanin.size() );

    node_num_fanin.push_back( num_fanin );

    for ( auto&& fo : fanouts )
    {
      add_fanin( fo, node_id );
    }

    if ( is_in_last_layer )
    {
      for ( auto slot = 0u; slot < num_fanin; slot++ )
      {
        last_layer_leaves.push_back( node_id );
      }
    }

    num_gates_of_fanin[num_fanin]++;

    return node_id;
  }

  size_t add_leaf_node( const part& fanouts = {} )
  {
    auto input_slot = add_internal_node( 0u, fanouts, false );
    input_slots.push_back( input_slot );
    return input_slot;
  }

  /*
   * Make a copy  with empty lists for last_layer_leaves, other_leaves, and input_slots.
   */
  aqfp_logical_network_t copy_without_leaves() const
  {
    aqfp_logical_network_t res;

    res.is_partial_dag = is_partial_dag;
    res.is_dag = is_dag;
    res.level = level;

    res.nodes = nodes;
    res.node_num_fanin = node_num_fanin;
    res.num_gates_of_fanin = num_gates_of_fanin;

    return res;
  }

  aqfp_logical_network_t copy_with_last_layer_leaves()
  {
    auto res = copy_without_leaves();
    res.last_layer_leaves = last_layer_leaves;
    return res;
  }

  /*
   * Create a aqfp_logical_structure with a single gate.
   */
  static aqfp_logical_network_t get_root( size_t num_fanin )
  {
    aqfp_logical_network_t net;

    net.is_partial_dag = true;
    net.is_dag = false;
    net.level = 1u;

    net.add_internal_node( num_fanin, {}, true );

    return net;
  }
};

using aqfp_logical_network = aqfp_logical_network_t<>;

std::string as_string( const aqfp_logical_network& net )
{
  std::vector<std::string> nodes;
  for ( auto i = 0u; i < net.nodes.size(); i++ )
  {
    std::vector<int> temp = net.nodes[i];
    while ( temp.size() < net.node_num_fanin[i] )
    {
      temp.push_back( i );
    }
    nodes.push_back( fmt::format( "{} = ({})", i, fmt::join( temp, " " ) ) );
  }
  return fmt::format( "{} N = {} ID {} DAG [ {} ] I = [{}] LL = [{}] OL = [{}] Z = {}", ( net.is_partial_dag ? "pdag" : ( net.is_dag ? "dag" : "invalid" ) ),
                      net.num_gates(),
                      net.pdag_id,
                      fmt::join( nodes, " " ),
                      fmt::join( net.input_slots, " " ),
                      fmt::join( net.last_layer_leaves, " " ),
                      fmt::join( net.other_leaves, " " ),
                      ( net.zero_input > 0 ) ? fmt::format( "{}", net.zero_input ) : "NA" );
}

/*
 * Compute the dag obtained from 'pdag' by combining slots according to 'p'.
 */
aqfp_logical_network get_next_dag( const aqfp_logical_network& orig, const partition& p )
{
  assert( orig.is_partial_dag );

  aqfp_logical_network net = orig.copy_without_leaves();

  net.is_partial_dag = false;
  net.is_dag = true;

  for ( auto&& q : p )
  {
    net.add_leaf_node( q );
  }

  return net;
}

/*
 * Compute the aqpf_logical_network obtained by combining the slots of 'orig' as
 * indicated by 'p'.
 * TODO: generate all pdags with different gate types
 */
aqfp_logical_network get_next_partial_dag( const aqfp_logical_network& orig, const partition& p,
                                           const std::vector<int>& other_leaves )
{
  assert( orig.is_partial_dag );

  aqfp_logical_network net = orig.copy_without_leaves();

  assert( net.is_partial_dag );
  assert( !net.is_dag );

  for ( auto&& q : p )
  {
    net.add_internal_node( 3u, q, true );
  }

  net.other_leaves = other_leaves; // set other leaves

  return net;
}

std::vector<aqfp_logical_network> add_node_recur( const aqfp_logical_network& orig, const std::vector<part>& p, size_t part_ind, std::map<size_t, size_t>& max_allowed_of_fanin, const dag_params& params )
{
  if ( part_ind == p.size() )
  {
    return { orig.copy_without_leaves() };
  }

  std::vector<aqfp_logical_network> res;

  // what fanin gate to use for part at part_ind?
  for ( auto&& fin : params.allowed_num_fanins )
  {
    if ( max_allowed_of_fanin[fin] == 0 )
    {
      continue;
    }
    max_allowed_of_fanin[fin]--;

    auto temp = add_node_recur( orig, p, part_ind + 1, max_allowed_of_fanin, params );

    for ( auto&& t : temp )
    {
      auto net = t.copy_with_last_layer_leaves();
      net.add_internal_node( fin, p[part_ind], true );
      res.push_back( net );
    }

    max_allowed_of_fanin[fin]++;
  }

  return res;
}

/*
 * Compute the aqpf_logical_network obtained by combining the slots of 'orig' as
 * indicated by 'p'.
 */
std::vector<aqfp_logical_network> get_next_partial_dags( const aqfp_logical_network& orig, const partition& p,
                                                         const std::vector<int>& other_leaves, const dag_params& params )
{

  auto max_allowed_of_fanin = params.max_gates_of_fanin;
  for ( auto it = orig.num_gates_of_fanin.begin(); it != orig.num_gates_of_fanin.end(); it++ )
  {
    assert( max_allowed_of_fanin[it->first] >= it->second );
    max_allowed_of_fanin[it->first] -= it->second;
  }

  std::vector<part> q( p.begin(), p.end() );

  auto res = add_node_recur( orig, q, 0, max_allowed_of_fanin, params );

  for ( auto&& net : res )
  {
    net.other_leaves = other_leaves;
  }

  return res;
}

std::vector<aqfp_logical_network> get_dags_from_partial_dag( const aqfp_logical_network& net, size_t leaves_lim = 0 )
{
  std::vector<int> leaves = net.last_layer_leaves;
  leaves.insert( leaves.end(), net.other_leaves.begin(), net.other_leaves.end() );

  std::sort( leaves.begin(), leaves.end() );

  auto max_counts = net.max_equal_fanins();

  auto partitions = detail::get_all_partitions( leaves, max_counts, leaves_lim, 0 /* unlimited part sizes (fanouts) */ );

  std::vector<aqfp_logical_network> result;
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

/*
 * extend the current aqfp logical network by one more level
 */
std::vector<aqfp_logical_network> get_layer_extension( const aqfp_logical_network& net, const dag_params& params )
{
  // choose which last layer slots to use
  // choose which non-last later slots to use
  // get partitions of the chosen last layer slots
  // extend the partitions using chosen non-last layer slots

  std::vector<aqfp_logical_network> result;

  auto max_counts = net.max_equal_fanins();

  auto last_options = detail::get_sub_lists( net.last_layer_leaves );
  auto other_options = detail::get_sub_lists( net.other_leaves );

  const std::map<int, size_t> last_counts = detail::get_frequencies( net.last_layer_leaves );
  const std::map<int, size_t> other_counts = detail::get_frequencies( net.other_leaves );

  for ( auto&& last : last_options )
  {
    if ( last.empty() )
      continue;

    // fmt::print( "last = [{}]\n", fmt::join( last, ", " ) );

    auto last_counts_cpy = last_counts;
    for ( auto&& e : last )
    {
      assert( last_counts_cpy[e] > 0 );
      last_counts_cpy[e]--;
    }

    for ( auto&& other : other_options )
    {
      // fmt::print( "\tother = [{}]\n", fmt::join( other, ", " ) );

      auto other_counts_cpy = other_counts;
      for ( auto&& e : other )
      {
        assert( other_counts_cpy[e] > 0 );
        other_counts_cpy[e]--;
      }

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
        auto last_layers_partitions = detail::get_all_partitions( last,
                                                                  max_counts,
                                                                  params.max_gates > 0 ? params.max_gates - net.num_gates() : 0u /* max parts */,
                                                                  params.max_num_fanout /* unlimited part size (fanouts) */ );
        for ( auto p : last_layers_partitions )
        {
          auto extensions = detail::extend_partitions( other,
                                                       p,
                                                       max_counts,
                                                       params.max_num_fanout /* unlimited part size (fanouts) */ );
          for ( auto q : extensions )
          {
            auto temp = get_next_partial_dags( net, q, other_leaves_new, params );
            for ( auto&& r : temp )
            {
              r.level++;
              result.push_back( r );
            }
          }
        }
      }
    }
  }

  return result;
}

template<typename Fn>
class gen_dag
{
public:
  gen_dag( const dag_params& params, Fn&& fn ) : params( params ), cost_fn( fn ),
                                                 pq( [&fn]( auto&& x, auto&& y ) -> bool { return fn( x ) > fn( y ) || ( fn( x ) == fn( y ) && y.is_dag && x.is_partial_dag ); } )
  {
    for ( auto&& fin : params.allowed_num_fanins )
    {
      if ( params.max_gates_of_fanin.at( fin ) > 0 )
      {
        auto root = aqfp_logical_network::get_root( fin );
        root.pdag_id = ( ++next_pdag_id );
        pq.push( root );
      }
    }
  }

  std::optional<aqfp_logical_network> next()
  {
    while ( true )
    {
      if ( pq.empty() )
      {
        return std::nullopt;
      }
      else
      {
        auto res = pq.top();
        pq.pop();

        if ( res.is_partial_dag )
        {
          if ( params.max_level > res.level )
          {
            auto ext = get_layer_extension( res, params );

            // std::for_each( std::execution::par_unseq, ext.begin(), ext.end(), []( auto&& item ) { cost_fn( item ); } );
            std::for_each( ext.begin(), ext.end(), [&]( auto&& item ) { cost_fn( item ); } );

            for ( auto&& t : ext )
            {
              t.pdag_id = ( ++next_pdag_id );
              pq.push( t );
            }
          }

          auto dags = get_dags_from_partial_dag( res, params.max_num_in );

          // std::for_each( std::execution::par_unseq, dags.begin(), dags.end(), []( auto&& items ) { cost_fn( item ); } );
          std::for_each( dags.begin(), dags.end(), [&]( auto&& item ) { cost_fn( item ); } );

          for ( auto&& t : dags )
          {
            t.pdag_id = res.pdag_id;
            pq.push( t );
          }
        }
        else
        {
          if ( res.input_slots.size() > params.max_num_in )
          {
            continue;
          }
          else
          {
            return res;
          }
        }
      }
    }
  }

  template<typename FeasibilityFn>
  std::optional<aqfp_logical_network> next_potentially_feasible( FeasibilityFn&& pdag_is_feasible )
  {
    while ( true )
    {
      if ( pq.empty() )
      {
        return std::nullopt;
      }
      else
      {
        auto res = pq.top();
        pq.pop();

        if ( res.is_partial_dag )
        {
          if ( params.max_level > res.level )
          {

            auto t0 = std::chrono::high_resolution_clock::now();

            auto ext = get_layer_extension( res, params );

            auto t1 = std::chrono::high_resolution_clock::now();

            // std::for_each( std::execution::par_unseq, ext.begin(), ext.end(), []( auto&& item ) { cost_fn( item ); } );
            auto block_size = 96u;

            for ( auto i = 0u; i < ext.size(); i += block_size )
            {
              std::vector<std::future<double>> temp;
              for ( auto j = i; j < i + block_size && j < ext.size(); j++ )
              {
                temp.push_back( std::async(
                    std::launch::async,
                    [&]( aqfp_logical_network& item ) { cost_fn( item ); return item.computed_cost; }, std::ref( ext[j] ) ) );
              }

              for ( auto j = i; j < i + block_size && j < ext.size(); j++ )
              {
                ext[j].computed_cost = temp[j - i].get();
              }
            }

            auto t2 = std::chrono::high_resolution_clock::now();

            auto d1 = std::chrono::duration_cast<std::chrono::microseconds>( t1 - t0 );
            auto d2 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );

            fmt::print( "generating {} partial dags took {} microseconds and computing costs took {} microseconds\n", ext.size(), d1.count(), d2.count() );

            for ( auto&& t : ext )
            {
              t.pdag_id = ( ++next_pdag_id );
              pq.push( t );
            }
          }
          if ( pdag_is_feasible( res ) )
          {
            auto t0 = std::chrono::high_resolution_clock::now();

            auto dags = get_dags_from_partial_dag( res, params.max_num_in );

            auto t1 = std::chrono::high_resolution_clock::now();

            // std::for_each( std::execution::par_unseq, dags.begin(), dags.end(), []( auto&& items ) { cost_fn( item ); } );
            auto block_size = 96u;

            for ( auto i = 0u; i < dags.size(); i += block_size )
            {
              std::vector<std::future<double>> temp;
              for ( auto j = i; j < i + block_size && j < dags.size(); j++ )
              {
                temp.push_back( std::async(
                    std::launch::async,
                    [&]( aqfp_logical_network& item ) { cost_fn( item ); return item.computed_cost; }, std::ref( dags[j] ) ) );
              }

              for ( auto j = i; j < i + block_size && j < dags.size(); j++ )
              {
                dags[j].computed_cost = temp[j - i].get();
              }
            }

            auto t2 = std::chrono::high_resolution_clock::now();

            auto d1 = std::chrono::duration_cast<std::chrono::microseconds>( t1 - t0 );
            auto d2 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 );

            fmt::print( "generating {} dags took {} microseconds and computing costs took {} microseconds\n", dags.size(), d1.count(), d2.count() );

            for ( auto&& t : dags )
            {
              t.pdag_id = res.pdag_id;
              pq.push( t );
            }
          }
        }
        else
        {
          if ( res.input_slots.size() > params.max_num_in )
          {
            continue;
          }
          else
          {
            return res;
          }
        }
      }
    }
  }

private:
  const dag_params params;
  std::function<double( aqfp_logical_network& )> cost_fn;
  std::priority_queue<aqfp_logical_network,
                      std::vector<aqfp_logical_network>, std::function<double( aqfp_logical_network&, aqfp_logical_network& )>>
      pq;
  size_t next_pdag_id = 0u;
};

} // namespace mockturtle
