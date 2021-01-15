#pragma once

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <vector>

#include <fmt/format.h>

namespace mockturtle
{

/**
 * Represents a single-output connected Partial DAG or DAG. A partial DAG is a network with majority gates
 * where some gates may have unconnected fanin slots.
 * A DAG is a network with majority gates obtained from a partial DAG by specifying
 * which unconnected fanin slots connect to the same primary input.
 * Optionally, a DAG may designate which fanin slots are connected the constant 0.
 * Unlike logic networks elsewhere in mockturtle, gates are numbered from 0 starting from the top gate.
 */
template<typename _NodeT = int>
struct aqfp_logical_network_t
{
  using NodeT = _NodeT;

  std::vector<uint32_t> node_num_fanin;  // number of fanins of each node
  std::vector<std::vector<NodeT>> nodes; // fanins of nodes
  std::vector<NodeT> input_slots;        // identifiers of the input slots (bundles of fanins where the inputs will be connected)
  NodeT zero_input = 0;                  // id of the input slot that is connected to contant 0

  /* meta data */
  uint32_t pdag_id = 0u;                           // an ID unique to each partial DAG
  bool is_partial_dag = false;                     // is this a partial DAG?
  uint32_t num_levels = 0u;                        // current number of levels
  double computed_cost = -1.0;                     // computed cost or -1.0 if not
  std::map<uint32_t, uint32_t> num_gates_of_fanin; // number of gates of each size
  std::vector<NodeT> last_layer_leaves;            // remaining fanin slots in the last layer
  std::vector<NodeT> other_leaves;                 // remaining fanin slots of the other layers

  /**
   * Compare to logical networks for equality.
   */
  bool operator==( const aqfp_logical_network_t<NodeT>& rhs ) const
  {
    if ( node_num_fanin != rhs.node_num_fanin )
      return false;

    if ( nodes.size() != rhs.nodes.size() )
      return false;
    for ( auto i = 0u; i < nodes.size(); i++ )
    {
      auto x1 = nodes[i];
      auto x2 = rhs.nodes[i];
      std::sort( x1.begin(), x1.end() );
      std::sort( x2.begin(), x2.end() );
      if ( x1 != x2 )
        return false;
    }

    auto y1 = input_slots;
    auto y2 = rhs.input_slots;
    std::sort( y1.begin(), y1.end() );
    std::sort( y2.begin(), y2.end() );
    if ( y1 != y2 )
      return false;

    if ( zero_input != rhs.zero_input )
      return false;

    return true;
  }

  /**
   * Decode a string representation of a DAG into a DAG.
   * Format: ng ni zi k0 g0f0 g0f1 .. g0fk0 k1 g1f0 g1f1 .. g1fk1 ....
   * ng := num gates, ni := num inputs, zi := zero input, ki := num fanin of i-th gate, gifj = j-th fanin of i-th gate
   */
  void decode_dag( std::string str )
  {
    pdag_id = 0;
    is_partial_dag = 0;
    num_levels = 0;
    computed_cost = -1.0;

    std::istringstream iss( str );
    auto ng = 0u;
    auto ni = 0u;
    auto zi = 0u;

    iss >> ng >> ni >> zi;
    zero_input = zi;

    std::vector<uint32_t> level( ng + ni, 0u );
    level[0u] = 1u;
    for ( auto i = 0u; i < ng; i++ )
    {
      auto nf = 0u;
      iss >> nf;

      node_num_fanin.push_back( nf );
      nodes.push_back( {} );

      num_gates_of_fanin[nf]++;

      for ( auto j = 0u; j < nf; j++ )
      {
        auto t = 0u;
        iss >> t;
        nodes[i].push_back( t );
        level[t] = std::max( level[t], level[i] + 1 );
      }
    }

    num_levels = *std::max_element( level.begin(), level.end() );

    for ( auto i = 0u; i < ni; i++ )
    {
      nodes.push_back( {} );
      node_num_fanin.push_back( 0 );
      input_slots.push_back( nodes.size() - 1 );
    }
  }

  /**
   * Encode a DAG as a string.
   * Format: ng ni zi k0 g0f0 g0f1 .. g0fk0 k1 g1f0 g1f1 .. g1fk1 ....
   * ng := num gates, ni := num inputs, zi := zero input, ki := num fanin of i-th gate, gifj = j-th fanin of i-th gate
   */
  std::string encode_as_string() const
  {
    assert( !is_partial_dag );

    std::stringstream ss;
    ss << num_gates() << " " << input_slots.size() << " " << zero_input;
    for ( auto i = 0u; i < num_gates(); i++ )
    {
      ss << fmt::format( " {} {}", node_num_fanin[i], fmt::join( nodes[i], " " ) );
    }

    return ss.str();
  }

  /**
   * Return the number of majority gates.
   */
  uint32_t num_gates() const
  {
    return nodes.size() - input_slots.size();
  }

  /**
   * Returns a map with maximum number of equal fanins that a gate may have without that gate being redundant.
   * A 3-input majority gate may not have any equal fanins as it would simplify otherwise.
   * A 5-input majoirty gate may have up to 2 equal fanins but if it had more, then it would simplify.
   */
  std::vector<uint32_t> max_equal_fanins() const
  {
    std::vector<uint32_t> res( num_gates() );
    for ( auto i = 0u; i < num_gates(); i++ )
    {
      res[i] = node_num_fanin[i] / 2;
    }
    return res;
  }

  /**
   * Add "fanin" as a fanin of "node".
   */
  void add_fanin( NodeT node, NodeT fanin )
  {
    nodes[node].push_back( fanin );
  }

  /**
   * Adds a new node with "num_fanin" fanins that is connected to "fanouts".
   * Optionally, it can be specified as a last layer node if at least one of its fanouts previously belonged to
   * the last layer.
   */
  uint32_t add_internal_node( uint32_t num_fanin = 3u, const std::multiset<NodeT>& fanouts = {}, bool is_in_last_layer = true )
  {
    uint32_t node_id = nodes.size();
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

    if ( num_fanin > 0u )
      num_gates_of_fanin[num_fanin]++;

    return node_id;
  }

  /**
   * Adds a new input node connected to "fanouts".
   */
  uint32_t add_leaf_node( const std::multiset<NodeT>& fanouts = {} )
  {
    /*
	! Change this to not call 'add_internal_nodes' and to not add empty nodes.
	! Might have to change everywhere we use nodes.size()
	*/
    auto input_slot = add_internal_node( 0u, fanouts, false );
    input_slots.push_back( input_slot );
    return input_slot;
  }

  /**
   * Make a copy  with empty last_layer_leaves and other_leaves.
   */
  aqfp_logical_network_t copy_without_leaves() const
  {
    aqfp_logical_network_t res{
        node_num_fanin,
        nodes,
        input_slots,
        zero_input,

        pdag_id,
        is_partial_dag,
        num_levels,
        -1.0,
        num_gates_of_fanin,
        {},
        {},
    };

    return res;
  }

  /**
   * Make a copy  with empty other_leaves.
   */
  aqfp_logical_network_t copy_with_last_layer_leaves()
  {
    aqfp_logical_network_t res{
        node_num_fanin,
        nodes,
        input_slots,
        zero_input,

        pdag_id,
        is_partial_dag,
        num_levels,
        -1.0,
        num_gates_of_fanin,
        last_layer_leaves,
        {},
    };

    return res;
  }

  /**
   * Create a aqfp_logical_structure with a single gate.
   */
  static aqfp_logical_network_t get_root( uint32_t num_fanin )
  {
    aqfp_logical_network_t net;

    net.is_partial_dag = true;
    net.num_levels = 1u;

    std::multiset<NodeT> fanouts = {};
    net.add_internal_node( num_fanin, fanouts, true );

    return net;
  }
};

/**
 * Returns a more readable string representation of a logical network.
 */
template<typename NodeT>
std::string as_string( const aqfp_logical_network_t<NodeT>& net )
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
  return fmt::format( "{} N = {} ID {} DAG [ {} ] I = [{}] LL = [{}] OL = [{}] Z = {}", ( net.is_partial_dag ? "PDAG" : "DAG" ),
                      net.num_gates(),
                      net.pdag_id,
                      fmt::join( nodes, " " ),
                      fmt::join( net.input_slots, " " ),
                      fmt::join( net.last_layer_leaves, " " ),
                      fmt::join( net.other_leaves, " " ),
                      ( net.zero_input > 0 ) ? fmt::format( "{}", net.zero_input ) : "NA" );
}

} // namespace mockturtle
