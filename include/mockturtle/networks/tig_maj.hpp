#pragma once

#include "tig.hpp"

namespace mockturtle
{

template <>
struct compute_function<three_input_function::majority> {
	template<typename T> T operator()( T a, T b, T c ) {
		return kitty::ternary_majority(a, b, c);
	}	   
};

using mig_signal = tig_network<three_input_function::majority>::signal;
using mig_network2 = tig_network<three_input_function::majority>;

template<> kitty::dynamic_truth_table mig_network2::node_function( const node& n ) const
{
    (void)n;
    kitty::dynamic_truth_table _maj( 3 );
    _maj._bits[0] = 0xe8;
    return _maj;
}

template<> bool mig_network2::is_maj( node const& n ) const
{
    return n > 0 && !is_ci( n );
}

template<> mig_signal mig_network2::create_maj( mig_signal a, mig_signal b, mig_signal c )
{
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

    /* trivial cases */
    if ( a.index == b.index )
    {
      return ( a.complement == b.complement ) ? a : c;
    }
    else if ( b.index == c.index )
    {
      return ( b.complement == c.complement ) ? b : a;
    }

    /*  complemented edges minimization */
    auto node_complement = false;
    if ( static_cast<unsigned>( a.complement ) + static_cast<unsigned>( b.complement ) +
             static_cast<unsigned>( c.complement ) >=
         2u )
    {
      node_complement = true;
      a.complement = !a.complement;
      b.complement = !b.complement;
      c.complement = !c.complement;
    }

    storage::element_type::node_type node;
    node.children[0] = a;
    node.children[1] = b;
    node.children[2] = c;

    /* structural hashing */
    const auto it = _storage->hash.find( node );
    if ( it != _storage->hash.end() )
    {
      return {it->second, node_complement};
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

    return {index, node_complement};
}

template<> mig_signal mig_network2::create_and( mig_signal a, mig_signal b)
{
	return create_maj( get_constant(false), a, b);
}

template<> mig_signal mig_network2::clone_node(mig_network2 const& other, node const& source, std::vector<mig_signal> const& children ) {
	(void)other;
	(void)source;
	assert( children.size() == 3u );
	return create_maj( children[0u], children[1u], children[2u] );
}
template<> mig_signal mig_network2::create_network_node( mig_signal a, mig_signal b, mig_signal c ) {
	return create_maj( a, b, c);
}

} // namespace mockturtle end
