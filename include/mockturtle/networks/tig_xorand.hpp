#pragma once

#include "tig.hpp"

namespace mockturtle
{


template <>
struct compute_function<three_input_function::xorand> {
	template<typename T> T operator()( T a, T b, T c ) {
    	return ternary_operation(a, b, c, [](auto a, auto b, auto c) { return a & (b ^ c); });
	}	   
};

using xaig_signal = tig_network<three_input_function::xorand>::signal;
using xaig_network = tig_network<three_input_function::xorand>;


template<> xaig_signal xaig_network::create_and( xaig_signal a, xaig_signal b);
template<> xaig_signal xaig_network::create_xor( xaig_signal a, xaig_signal b);


template<> kitty::dynamic_truth_table xaig_network::node_function( const node& n ) const
{
    (void)n;
    kitty::dynamic_truth_table tt( 3 );
    tt._bits[0] = 0x28;
    return tt;
}

template<> bool xaig_network::is_xorand( node const& n ) const
{
	return n > 0 && !is_ci( n );
}

template<> xaig_signal xaig_network::create_xorand( xaig_signal a, xaig_signal b, xaig_signal c )
{		
	/* xorand(a b c) = a (b xor z)*/

	/* if some inputs are constants */
	if (a.index == 0) { 					// (0 b c) = 0, (1 b c) = b xor c
		return (a.complement == 0) ? get_constant( false ) : create_xor(b, c); 
	} else {
		if (b.index == 0) { 				// (a 0 c) = a c, (a 1 c) = a c'
			return (b.complement == 0) ? create_and(a, c) : create_and(a, !c);
		} else {
			if (c.index == 0) { 			// (a b 0) = a b,  (a b 1) = a b'
				return (c.complement == 0) ? create_and(a, b) : create_and(a, !b);
			}
		}
	}

	/* if there are repeated inputs */
	if (a.index == b.index) {   			// (a a c) = a c', (a a' c) = a c
		return (a.complement == b.complement) ? create_and(a, !c) : create_and(a, c);
	} else {
		if (a.index == c.index) {			// (a b a) = a b', (a b a') = a b 
			return (a.complement == c.complement) ? create_and(a, !b) : create_and(a, b);
		} else {
			if (b.index == c.index) {		// (a b b) = 0, (a b b') = a
				return (b.complement == c.complement) ? get_constant( false ) : a;
			}
		}
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

template<> xaig_signal xaig_network::create_and( xaig_signal a, xaig_signal b)
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
    node.children[0] = a;
    node.children[1] = b;
    node.children[2] = get_constant( false );

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


template<> xaig_signal xaig_network::create_xor( xaig_signal a, xaig_signal b)
{
	if (a.index > b.index) {
		std::swap(a, b);
	}

	if (a.index == 0) { 						// (0 b) = b, (1 b) = b'
		return (a.complement == 0) ? b : !b;	
	} 

	if (a.index == b.index) {
		return (a.complement == b.complement) ? get_constant( false ) : !get_constant( false );
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



template<> xaig_signal xaig_network::clone_node( xaig_network const& other, node const& source, std::vector<xaig_signal> const& children ) {
	(void)other;
	(void)source;
	assert( children.size() == 3u );
	return create_xorand( children[0u], children[1u], children[2u] );
}


template<> xaig_signal xaig_network::create_network_node( xaig_signal a, xaig_signal b, xaig_signal c ) {
	return create_xorand( a, b, c);
}


} // namespace mockturtle end
