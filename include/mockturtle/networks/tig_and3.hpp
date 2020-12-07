#pragma once

#include "tig.hpp"

namespace mockturtle
{

template <>
struct compute_function<three_input_function::and3> {
	template<typename T> T operator()( T a, T b, T c ) {
		return ternary_operation(a, b, c, [](auto a, auto b, auto c) {return (a & b & c);});
	}	   
};

using a3ig_signal = tig_network<three_input_function::and3>::signal;
using a3ig_network = tig_network<three_input_function::and3>;

template<> a3ig_signal a3ig_network::create_and( a3ig_signal a, a3ig_signal b);

template<> kitty::dynamic_truth_table a3ig_network::node_function( const node& n ) const
{
    (void)n;
    kitty::dynamic_truth_table tt( 3 );
    tt._bits[0] = 0x80;
    return tt;
}

template<> bool a3ig_network::is_and3( node const& n ) const
{
    return n > 0 && !is_ci( n );
}

template<> a3ig_signal a3ig_network::create_and3( a3ig_signal a, a3ig_signal b, a3ig_signal c )
{
	
	/* and3(a b c) = a b c */
    /* order inputs */
    if ( a.index > b.index )
    {
      std::swap( a, b );
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }
    else
    {
      if ( b.index > c.index )
        std::swap( b, c );
      if ( a.index > b.index )
        std::swap( a, b );
    }

	/* if some inputs are constants */
	if (a.index == 0) { 					// (0 b c) = 0, (1 b c) = b c
		return (a.complement == 0) ? get_constant( false ) : create_and(b, c); 
	} // If a is not constant b or c cannot be constants (sorted inputs)

	/* if there are repeated inputs */
	if (a.index == b.index) {   			// (a a c) = a c, (a a' c) = 0
		return (a.complement == b.complement) ? create_and(a, c) : get_constant( false );
	} 
	
	if (b.index == c.index) {		// (a b b) = a b, (a b b') = 0
		return (b.complement == c.complement) ? create_and(a, b) : get_constant( false );
	}

    storage::element_type::node_type node;
    node.children[0] = a;
    node.children[1] = b;
    node.children[2] = c;

    /* structural hashing */
    const auto it = _storage->hash.find( node );
    if ( it != _storage->hash.end() )
    {
      return {it->second, 0};
    }

    const auto index = _storage->nodes.size();

    if ( index >= .9 * _storage->nodes.capacity() )
    {
      _storage->nodes.reserve( static_cast<uint64_t>( 3.1415f * index ) );
      _storage->hash.reserve( static_cast<uint64_t>( 3.1415f * index ) );
    }

    _storage->nodes.push_back( node );

    _storage->hash[node] = index;

    /* increase ref-count to children */
    _storage->nodes[a.index].data[0].h1++;
    _storage->nodes[b.index].data[0].h1++;
    _storage->nodes[c.index].data[0].h1++;

    for ( auto const& fn : _events->on_add )
    {
      fn( index );
    }

    return {index, 0};
}

template<> a3ig_signal a3ig_network::create_and( a3ig_signal a, a3ig_signal b)
{
	if (a.index > b.index) {
		std::swap(a, b);
	}

	if (a.index == 0) { 			// (0 b) = 0, (1 b) = b
		return (a.complement == 0) ? get_constant( false ) : b;	
	} 

	if (a.index == b.index) {
		return (a.complement == b.complement) ? a : get_constant( false );
	}


    storage::element_type::node_type node;
    node.children[0] = !get_constant( false );
    node.children[1] = a;
    node.children[2] = b;

    /* structural hashing */
    const auto it = _storage->hash.find( node );
    if ( it != _storage->hash.end() )
    {
      return {it->second, 0};
    }

    const auto index = _storage->nodes.size();

    if ( index >= .9 * _storage->nodes.capacity() )
    {
      _storage->nodes.reserve( static_cast<uint64_t>( 3.1415f * index ) );
      _storage->hash.reserve( static_cast<uint64_t>( 3.1415f * index ) );
    }

    _storage->nodes.push_back( node );

    _storage->hash[node] = index;

    /* increase ref-count to children */
    _storage->nodes[a.index].data[0].h1++;
    _storage->nodes[b.index].data[0].h1++;
    _storage->nodes[0].data[0].h1++;

    for ( auto const& fn : _events->on_add )
    {
      fn( index );
    }

    return {index, 0};
}

template<> a3ig_signal a3ig_network::clone_node(a3ig_network const& other, node const& source, std::vector<a3ig_signal> const& children ) {
	(void)other;
	(void)source;
	assert( children.size() == 3u );
	return create_and3( children[0u], children[1u], children[2u] );
}

template<> a3ig_signal a3ig_network::create_network_node( a3ig_signal a, a3ig_signal b, a3ig_signal c ) {
	return create_and3( a, b, c);
}

} // namespace mockturtle end
