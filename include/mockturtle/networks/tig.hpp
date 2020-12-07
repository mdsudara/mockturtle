/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file tig.hpp
  \brief MIG logic network implementation

  \author Eleonora Testa, Dewmini Sudara
  */

#pragma once

#include <memory>
#include <optional>
#include <stack>
#include <string>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operators.hpp>

#include "../traits.hpp"
#include "../utils/algorithm.hpp"
#include "detail/foreach.hpp"
#include "events.hpp"
#include "storage.hpp"

namespace mockturtle
{
	enum class three_input_function { dot, onehot, mux, andxor, xorand, gamble, orand, majority, and3, xor3 };

	template <three_input_function Fn>
	struct compute_function {
		template<typename T>  T operator()( T a, T b, T c ); 
	};

	struct tig_storage_data {
		uint32_t num_pis = 0u;
		uint32_t num_pos = 0u;
		std::vector<int8_t> latches;
		uint32_t trav_id = 0u;
	};

	/*! \brief T-Inverter-Graph storage container

	  TIGs have nodes with fan-in 3.  We split of one bit of the index pointer to
	  store a complemented attribute.  Every node has 64-bit of additional data
	  used for the following purposes:

	  `data[0].h1`: Fan-out size (we use MSB to indicate whether a node is dead)
	  `data[0].h2`: Application-specific value
	  `data[1].h1`: Visited flag
	  */

	using three_input_func =  bool (*)(bool, bool, bool);

	using tig_node = regular_node<3, 2, 1>;
	using tig_storage = storage<tig_node, tig_storage_data>;

	struct tig_signal
	{
		tig_signal() = default;

		tig_signal( uint64_t index, uint64_t complement )
			: complement( complement ), index( index )
		{
		}

		explicit tig_signal( uint64_t data )
			: data( data )
		{
		}

		tig_signal( tig_storage::node_type::pointer_type const& p )
			: complement( p.weight ), index( p.index )
		{
		}

		union {
			struct
			{
				uint64_t complement : 1;
				uint64_t index : 63;
			};
			uint64_t data;
		};

		tig_signal operator!() const
		{
			return tig_signal( data ^ 1 );
		}

		tig_signal operator+() const
		{
			return {index, 0};
		}

		tig_signal operator-() const
		{
			return {index, 1};
		}

		tig_signal operator^( bool complement ) const
		{
			return tig_signal( data ^ ( complement ? 1 : 0 ) );
		}

		bool operator==( tig_signal const& other ) const
		{
			return data == other.data;
		}

		bool operator!=( tig_signal const& other ) const
		{
			return data != other.data;
		}

		bool operator<( tig_signal const& other ) const
		{
			return data < other.data;
		}

		operator tig_storage::node_type::pointer_type() const
		{
			return {index, complement};
		}
	};


	template <three_input_function NodeFunc>
		class tig_network
		{
			public:
#pragma region Types and constructors
				static constexpr auto min_fanin_size = 3u;
				static constexpr auto max_fanin_size = 3u;

				using base_type = tig_network<NodeFunc>;
				using storage = std::shared_ptr<tig_storage>;
				using node = uint64_t;
				using signal = tig_signal;

				tig_network()
					: _storage( std::make_shared<tig_storage>() ),
					_events( std::make_shared<typename decltype( _events )::element_type>() )
			{
			}

				tig_network( std::shared_ptr<tig_storage> storage )
					: _storage( storage ),
					_events( std::make_shared<typename decltype( _events )::element_type>() )
			{
			}
#pragma endregion

#pragma region Primary I / O and constants
				signal get_constant( bool value ) const
				{
					return {0, static_cast<uint64_t>( value ? 1 : 0 )};
				}

				signal create_pi( std::string const& name = std::string() )
				{
					(void)name;

					const auto index = _storage->nodes.size();
					auto& node = _storage->nodes.emplace_back();
					node.children[0].data = node.children[1].data = node.children[2].data = ~static_cast<uint64_t>( 0 );
					_storage->inputs.emplace_back( index );
					++_storage->data.num_pis;
					return {index, 0};
				}

				uint32_t create_po( signal const& f, std::string const& name = std::string() )
				{
					(void)name;

					/* increase ref-count to children */
					_storage->nodes[f.index].data[0].h1++;
					auto const po_index = static_cast<uint32_t>( _storage->outputs.size() );
					_storage->outputs.emplace_back( f.index, f.complement );
					++_storage->data.num_pos;
					return po_index;
				}

				signal create_ro( std::string const& name = std::string() )
				{
					(void)name;

					auto const index = _storage->nodes.size();
					auto& node = _storage->nodes.emplace_back();
					node.children[0].data = node.children[1].data = node.children[2].data = _storage->inputs.size();
					_storage->inputs.emplace_back( index );
					return {index, 0};
				}

				uint32_t create_ri( signal const& f, int8_t reset = 0, std::string const& name = std::string() )
				{
					(void)name;

					/* increase ref-count to children */
					_storage->nodes[f.index].data[0].h1++;
					auto const ri_index = static_cast<uint32_t>( _storage->outputs.size() );
					_storage->outputs.emplace_back( f.index, f.complement );
					_storage->data.latches.emplace_back( reset );
					return ri_index;
				}

				int8_t latch_reset( uint32_t index ) const
				{
					assert( index < _storage->data.latches.size() );
					return _storage->data.latches[index];
				}

				bool is_combinational() const
				{
					return ( static_cast<uint32_t>( _storage->inputs.size() ) == _storage->data.num_pis &&
							static_cast<uint32_t>( _storage->outputs.size() ) == _storage->data.num_pos );
				}

				bool is_constant( node const& n ) const
				{
					return n == 0;
				}

				bool is_ci( node const& n ) const
				{
					return _storage->nodes[n].children[0].data == _storage->nodes[n].children[1].data && _storage->nodes[n].children[0].data == _storage->nodes[n].children[2].data;
				}

				bool is_pi( node const& n ) const
				{
					return _storage->nodes[n].children[0].data == ~static_cast<uint64_t>( 0 ) && _storage->nodes[n].children[1].data == ~static_cast<uint64_t>( 0 ) && _storage->nodes[n].children[2].data == ~static_cast<uint64_t>( 0 );
				}

				bool is_ro( node const& n ) const
				{
					return _storage->nodes[n].children[0].data == _storage->nodes[n].children[1].data && _storage->nodes[n].children[0].data == _storage->nodes[n].children[2].data && _storage->nodes[n].children[0].data >= static_cast<uint64_t>( _storage->data.num_pis );
				}

				bool constant_value( node const& n ) const
				{
					(void)n;
					return false;
				}
#pragma endregion

#pragma region Create unary functions
				signal create_buf( signal const& a )
				{
					return a;
				}

				signal create_not( signal const& a )
				{
					return !a;
				}
#pragma endregion

#pragma region Create binary / ternary functions

				/* Remark: must be specialized. */
				signal create_and( signal a, signal b );

				signal create_nand( signal a, signal b )
				{
					return !create_and( a, b );
				}

				signal create_or( signal a, signal b )
				{
					return !create_and( !a, !b );
				}

				signal create_nor( signal a, signal b )
				{
					return !create_or( a, b );
				}

				signal create_lt( signal a, signal b )
				{
					return create_and( !a, b );
				}

				signal create_le( signal a, signal b )
				{
					return !create_and( a, !b );
				}

				signal create_network_node( signal a, signal b, signal c );

				signal create_dot( signal a, signal b, signal c);

				signal create_onehot( signal a, signal b, signal c);
				
				signal create_and3( signal a, signal b, signal c);

				signal create_gamble( signal a, signal b, signal c);
				
				signal create_andxor( signal a, signal b, signal c);
				
				signal create_xorand( signal a, signal b, signal c);
				
				signal create_orand( signal a, signal b, signal c);
				
				signal create_mux( signal a, signal b, signal c);

				/* Remark: to be specialized. */
				signal create_xor( signal a, signal b )
				{
					const auto fcompl = a.complement ^ b.complement;
					const auto c1 = create_and( +a, -b );
					const auto c2 = create_and( +b, -a );
					return create_and( !c1, !c2 ) ^ !fcompl;
				}

				/* Remark: to be specialized. */
				signal create_maj( signal a, signal b, signal c )
				{
					return create_or( create_and( a, b), create_and( c, create_or( a, b ) ) );
				} 

				signal create_ite( signal cond, signal f_then, signal f_else )
				{
					bool f_compl{false};
					if ( f_then.index < f_else.index )
					{
						std::swap( f_then, f_else );
						cond.complement ^= 1;
					}
					if ( f_then.complement )
					{
						f_then.complement = 0;
						f_else.complement ^= 1;
						f_compl = true;
					}

					return create_and( !create_and( !cond, f_else ), !create_and( cond, f_then ) ) ^ !f_compl;
				}

				signal create_xor3( signal const& a, signal const& b, signal const& c )
				{
					return create_xor( a, create_xor( b, c));
				}
#pragma endregion

#pragma region Create nary functions
				signal create_nary_and( std::vector<signal> const& fs )
				{
					return tree_reduce( fs.begin(), fs.end(), get_constant( true ), [this]( auto const& a, auto const& b ) { return create_and( a, b ); } );
				}

				signal create_nary_or( std::vector<signal> const& fs )
				{
					return tree_reduce( fs.begin(), fs.end(), get_constant( false ), [this]( auto const& a, auto const& b ) { return create_or( a, b ); } );
				}

				signal create_nary_xor( std::vector<signal> const& fs )
				{
					return tree_reduce( fs.begin(), fs.end(), get_constant( false ), [this]( auto const& a, auto const& b ) { return create_xor( a, b ); } );
				}
#pragma endregion

#pragma region Create arbitrary functions
				/* Remark: must be specialized. */
				signal clone_node( tig_network const& other, node const& source, std::vector<signal> const& children );
#pragma endregion

#pragma region Restructuring

				/* Remark: to be specialized. */
				std::optional<std::pair<node, signal>> replace_in_node( node const& n, node const& old_node, signal new_signal )
				{
					auto& node = _storage->nodes[n];

					signal child[3];

					bool found = false;
					for ( auto i = 0u; i < 3u; ++i )
					{
						if ( node.children[i].index == old_node )
						{
							found = true;
							new_signal.complement ^= node.children[i].weight;
							child[i] = new_signal;
						} else {
							child[i] = node.children[i];
						}
					}

					if ( !found )
					{
						return std::nullopt;
					}

					// node already in hash table
					storage::element_type::node_type _hash_obj;
					_hash_obj.children[0] = child[0];
					_hash_obj.children[1] = child[1];
					_hash_obj.children[2] = child[2];
					if ( const auto it = _storage->hash.find( _hash_obj ); it != _storage->hash.end() )
					{
						return std::make_pair( n, signal( it->second, 0 ) );
					}

					// remember before
					const auto old_child0 = signal{node.children[0]};
					const auto old_child1 = signal{node.children[1]};
					const auto old_child2 = signal{node.children[2]};

					// erase old node in hash table
					_storage->hash.erase( node );

					// insert updated node into hash table
					node.children[0] = child[0];
					node.children[1] = child[1];
					node.children[2] = child[2];
					_storage->hash[node] = n;

					// update the reference counter of the new signal
					_storage->nodes[new_signal.index].data[0].h1++;

					for ( auto const& fn : _events->on_modified )
					{
						fn( n, {old_child0, old_child1, old_child2} );
					}

					return std::nullopt;
				}

				void replace_in_outputs( node const& old_node, signal const& new_signal )
				{
					for ( auto& output : _storage->outputs )
					{
						if ( output.index == old_node )
						{
							output.index = new_signal.index;
							output.weight ^= new_signal.complement;

							// increment fan-in of new node
							_storage->nodes[new_signal.index].data[0].h1++;
						}
					}
				}

				void take_out_node( node const& n )
				{
					/* we cannot delete CIs or constants */
					if ( n == 0 || is_ci( n ) )
						return;

					auto& nobj = _storage->nodes[n];
					nobj.data[0].h1 = UINT32_C( 0x80000000 ); /* fanout size 0, but dead */
					_storage->hash.erase( nobj );

					for ( auto const& fn : _events->on_delete )
					{
						fn( n );
					}

					for ( auto i = 0u; i < 3u; ++i )
					{
						if ( fanout_size( nobj.children[i].index ) == 0 )
						{
							continue;
						}
						if ( decr_fanout_size( nobj.children[i].index ) == 0 )
						{
							take_out_node( nobj.children[i].index );
						}
					}
				}

				inline bool is_dead( node const& n ) const
				{
					return ( _storage->nodes[n].data[0].h1 >> 31 ) & 1;
				}

				void substitute_node( node const& old_node, signal const& new_signal )
				{
					std::stack<std::pair<node, signal>> to_substitute;
					to_substitute.push( {old_node, new_signal} );

					while ( !to_substitute.empty() )
					{
						const auto [_old, _new] = to_substitute.top();
						to_substitute.pop();

						for ( auto idx = 1u; idx < _storage->nodes.size(); ++idx )
						{
							if ( is_ci( idx ) )
								continue; /* ignore CIs */

							if ( const auto repl = replace_in_node( idx, _old, _new ); repl )
							{
								to_substitute.push( *repl );
							}
						}

						/* check outputs */
						replace_in_outputs( _old, _new );

						// reset fan-in of old node
						take_out_node( _old );
					}
				}

				void substitute_node_of_parents( std::vector<node> const& parents, node const& old_node, signal const& new_signal )
				{
					for ( auto& p : parents )
					{
						auto& n = _storage->nodes[p];
						for ( auto& child : n.children )
						{
							if ( child.index == old_node )
							{
								child.index = new_signal.index;
								child.weight ^= new_signal.complement;

								// increment fan-in of new node
								_storage->nodes[new_signal.index].data[0].h1++;

								// decrement fan-in of old node
								_storage->nodes[old_node].data[0].h1--;
							}
						}
					}

					/* check outputs */
					for ( auto& output : _storage->outputs )
					{
						if ( output.index == old_node )
						{
							output.index = new_signal.index;
							output.weight ^= new_signal.complement;

							// increment fan-in of new node
							_storage->nodes[new_signal.index].data[0].h1++;

							// decrement fan-in of old node
							_storage->nodes[old_node].data[0].h1--;
						}
					}
				}
#pragma endregion

#pragma region Structural properties
				auto size() const
				{
					return static_cast<uint32_t>( _storage->nodes.size() );
				}

				auto num_cis() const
				{
					return static_cast<uint32_t>( _storage->inputs.size() );
				}

				auto num_cos() const
				{
					return static_cast<uint32_t>( _storage->outputs.size() );
				}

				uint32_t num_latches() const
				{
					return static_cast<uint32_t>( _storage->data.latches.size() );
				}

				auto num_pis() const
				{
					return _storage->data.num_pis;
				}

				auto num_pos() const
				{
					return _storage->data.num_pos;
				}

				auto num_registers() const
				{
					assert( static_cast<uint32_t>( _storage->inputs.size() - _storage->data.num_pis ) == static_cast<uint32_t>( _storage->outputs.size() - _storage->data.num_pos ) );
					return static_cast<uint32_t>( _storage->inputs.size() - _storage->data.num_pis );
				}

				auto num_gates() const
				{
					return static_cast<uint32_t>( _storage->hash.size() );
				}

				uint32_t fanin_size( node const& n ) const
				{
					if ( is_constant( n ) || is_ci( n ) )
						return 0;
					return 3;
				}

				uint32_t fanout_size( node const& n ) const
				{
					return _storage->nodes[n].data[0].h1 & UINT32_C( 0x7FFFFFFF );
				}

				uint32_t incr_fanout_size( node const& n ) const
				{
					return _storage->nodes[n].data[0].h1++ & UINT32_C( 0x7FFFFFFF );
				}

				uint32_t decr_fanout_size( node const& n ) const
				{
					return --_storage->nodes[n].data[0].h1 & UINT32_C( 0x7FFFFFFF );
				}

				bool is_and( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_or( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_xor( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_maj( node const& n ) const
				{
					(void)n;
					return false;
				}

			
				bool is_mux( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_xor3( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_and3( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_dot( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_onehot( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_orand( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_xorand( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_andxor( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_gamble( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_ite( node const& n ) const
				{
					return is_mux( n );
				}


				bool is_nary_and( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_nary_or( node const& n ) const
				{
					(void)n;
					return false;
				}

				bool is_nary_xor( node const& n ) const
				{
					(void)n;
					return false;
				}
#pragma endregion

#pragma region Functional properties
				/* Remark: must be specialized */
				kitty::dynamic_truth_table node_function( const node& n ) const;
#pragma endregion

#pragma region Nodes and signals
				node get_node( signal const& f ) const
				{
					return f.index;
				}

				signal make_signal( node const& n ) const
				{
					return signal( n, 0 );
				}

				bool is_complemented( signal const& f ) const
				{
					return f.complement;
				}

				uint32_t node_to_index( node const& n ) const
				{
					return static_cast<uint32_t>( n );
				}

				node index_to_node( uint32_t index ) const
				{
					return index;
				}

				node ci_at( uint32_t index ) const
				{
					assert( index < _storage->inputs.size() );
					return *( _storage->inputs.begin() + index );
				}

				signal co_at( uint32_t index ) const
				{
					assert( index < _storage->outputs.size() );
					return *( _storage->outputs.begin() + index );
				}

				node pi_at( uint32_t index ) const
				{
					assert( index < _storage->data.num_pis );
					return *( _storage->inputs.begin() + index );
				}

				signal po_at( uint32_t index ) const
				{
					assert( index < _storage->data.num_pos );
					return *( _storage->outputs.begin() + index );
				}

				node ro_at( uint32_t index ) const
				{
					assert( index < _storage->inputs.size() - _storage->data.num_pis );
					return *( _storage->inputs.begin() + _storage->data.num_pis + index );
				}

				signal ri_at( uint32_t index ) const
				{
					assert( index < _storage->outputs.size() - _storage->data.num_pos );
					return *( _storage->outputs.begin() + _storage->data.num_pos + index );
				}

				uint32_t ci_index( node const& n ) const
				{
					assert( _storage->nodes[n].children[0].data == _storage->nodes[n].children[1].data &&
							_storage->nodes[n].children[0].data == _storage->nodes[n].children[2].data );
					return static_cast<uint32_t>( _storage->nodes[n].children[0].data );
				}

				uint32_t co_index( signal const& s ) const
				{
					uint32_t i = -1;
					foreach_co( [&]( const auto& x, auto index ) {
							if ( x == s )
							{
							i = index;
							return false;
							}
							return true;
							} );
					return i;
				}

				uint32_t pi_index( node const& n ) const
				{
					assert( _storage->nodes[n].children[0].data == _storage->nodes[n].children[1].data && 
							_storage->nodes[n].children[0].data == _storage->nodes[n].children[2].data);
					assert( _storage->nodes[n].children[0].data < _storage->data.num_pis );

					return static_cast<uint32_t>( _storage->nodes[n].children[0].data );
				}

				uint32_t po_index( signal const& s ) const
				{
					uint32_t i = -1;
					foreach_po( [&]( const auto& x, auto index ) {
							if ( x == s )
							{
							i = index;
							return false;
							}
							return true;
							} );
					return i;
				}

				uint32_t ro_index( node const& n ) const
				{
					assert( _storage->nodes[n].children[0].data == _storage->nodes[n].children[1].data && 
							_storage->nodes[n].children[0].data == _storage->nodes[n].children[2].data );
					assert( _storage->nodes[n].children[0].data >= _storage->data.num_pis );

					return static_cast<uint32_t>( _storage->nodes[n].children[0].data - _storage->data.num_pis );
				}

				uint32_t ri_index( signal const& s ) const
				{
					uint32_t i = -1;
					foreach_ri( [&]( const auto& x, auto index ) {
							if ( x == s )
							{
							i = index;
							return false;
							}
							return true;
							} );
					return i;
				}

				signal ro_to_ri( signal const& s ) const
				{
					return *( _storage->outputs.begin() + _storage->data.num_pos + _storage->nodes[s.index].children[0].data - _storage->data.num_pis );
				}

				node ri_to_ro( signal const& s ) const
				{
					return *( _storage->inputs.begin() + _storage->data.num_pis + ri_index( s ) );
				}
#pragma endregion

#pragma region Node and signal iterators
				template<typename Fn>
					void foreach_node( Fn&& fn ) const
					{
						auto r = range<uint64_t>( _storage->nodes.size() );
						detail::foreach_element_if(
								r.begin(), r.end(),
								[this]( auto n ) { return !is_dead( n ); },
								fn );
					}

				template<typename Fn>
					void foreach_ci( Fn&& fn ) const
					{
						detail::foreach_element( _storage->inputs.begin(), _storage->inputs.end(), fn );
					}

				template<typename Fn>
					void foreach_co( Fn&& fn ) const
					{
						detail::foreach_element( _storage->outputs.begin(), _storage->outputs.end(), fn );
					}

				template<typename Fn>
					void foreach_pi( Fn&& fn ) const
					{
						detail::foreach_element( _storage->inputs.begin(), _storage->inputs.begin() + _storage->data.num_pis, fn );
					}

				template<typename Fn>
					void foreach_po( Fn&& fn ) const
					{
						detail::foreach_element( _storage->outputs.begin(), _storage->outputs.begin() + _storage->data.num_pos, fn );
					}

				template<typename Fn>
					void foreach_ro( Fn&& fn ) const
					{
						detail::foreach_element( _storage->inputs.begin() + _storage->data.num_pis, _storage->inputs.end(), fn );
					}

				template<typename Fn>
					void foreach_ri( Fn&& fn ) const
					{
						detail::foreach_element( _storage->outputs.begin() + _storage->data.num_pos, _storage->outputs.end(), fn );
					}

				template<typename Fn>
					void foreach_register( Fn&& fn ) const
					{
						static_assert( detail::is_callable_with_index_v<Fn, std::pair<signal, node>, void> ||
								detail::is_callable_without_index_v<Fn, std::pair<signal, node>, void> ||
								detail::is_callable_with_index_v<Fn, std::pair<signal, node>, bool> ||
								detail::is_callable_without_index_v<Fn, std::pair<signal, node>, bool> );

						assert( _storage->inputs.size() - _storage->data.num_pis == _storage->outputs.size() - _storage->data.num_pos );
						auto ro = _storage->inputs.begin() + _storage->data.num_pis;
						auto ri = _storage->outputs.begin() + _storage->data.num_pos;
						if constexpr ( detail::is_callable_without_index_v<Fn, std::pair<signal, node>, bool> )
						{
							while ( ro != _storage->inputs.end() && ri != _storage->outputs.end() )
							{
								if ( !fn( std::make_pair( ri++, ro++ ) ) )
									return;
							}
						}
						else if constexpr ( detail::is_callable_with_index_v<Fn, std::pair<signal, node>, bool> )
						{
							uint32_t index{0};
							while ( ro != _storage->inputs.end() && ri != _storage->outputs.end() )
							{
								if ( !fn( std::make_pair( ri++, ro++ ), index++ ) )
									return;
							}
						}
						else if constexpr ( detail::is_callable_without_index_v<Fn, std::pair<signal, node>, void> )
						{
							while ( ro != _storage->inputs.end() && ri != _storage->outputs.end() )
							{
								fn( std::make_pair( *ri++, *ro++ ) );
							}
						}
						else if constexpr ( detail::is_callable_with_index_v<Fn, std::pair<signal, node>, void> )
						{
							uint32_t index{0};
							while ( ro != _storage->inputs.end() && ri != _storage->outputs.end() )
							{
								fn( std::make_pair( *ri++, *ro++ ), index++ );
							}
						}
					}

				template<typename Fn>
					void foreach_gate( Fn&& fn ) const
					{
						auto r = range<uint64_t>( 1u, _storage->nodes.size() ); // start from 1 to avoid constant
						detail::foreach_element_if(
								r.begin(), r.end(),
								[this]( auto n ) { return !is_ci( n ) && !is_dead( n ); },
								fn );
					}

				template<typename Fn>
					void foreach_fanin( node const& n, Fn&& fn ) const
					{
						if ( n == 0 || is_ci( n ) )
							return;

						static_assert( detail::is_callable_without_index_v<Fn, signal, bool> ||
								detail::is_callable_with_index_v<Fn, signal, bool> ||
								detail::is_callable_without_index_v<Fn, signal, void> ||
								detail::is_callable_with_index_v<Fn, signal, void> );

						// we don't use foreach_element here to have better performance
						if constexpr ( detail::is_callable_without_index_v<Fn, signal, bool> )
						{
							if ( !fn( signal{_storage->nodes[n].children[0]} ) )
								return;
							if ( !fn( signal{_storage->nodes[n].children[1]} ) )
								return;
							fn( signal{_storage->nodes[n].children[2]} );
						}
						else if constexpr ( detail::is_callable_with_index_v<Fn, signal, bool> )
						{
							if ( !fn( signal{_storage->nodes[n].children[0]}, 0 ) )
								return;
							if ( !fn( signal{_storage->nodes[n].children[1]}, 1 ) )
								return;
							fn( signal{_storage->nodes[n].children[2]}, 2 );
						}
						else if constexpr ( detail::is_callable_without_index_v<Fn, signal, void> )
						{
							fn( signal{_storage->nodes[n].children[0]} );
							fn( signal{_storage->nodes[n].children[1]} );
							fn( signal{_storage->nodes[n].children[2]} );
						}
						else if constexpr ( detail::is_callable_with_index_v<Fn, signal, void> )
						{
							fn( signal{_storage->nodes[n].children[0]}, 0 );
							fn( signal{_storage->nodes[n].children[1]}, 1 );
							fn( signal{_storage->nodes[n].children[2]}, 2 );
						}
					}
#pragma endregion

#pragma region Value simulation
template<typename Iterator>
	iterates_over_t<Iterator, bool>
	compute( node const& n, Iterator begin, Iterator end ) const
	{
		(void)end;

		assert( n != 0 && !is_ci( n ) );

		auto const& c1 = _storage->nodes[n].children[0];
		auto const& c2 = _storage->nodes[n].children[1];
		auto const& c3 = _storage->nodes[n].children[2];

		auto v1 = *begin++;
		auto v2 = *begin++;
		auto v3 = *begin++;

		return compute_function<NodeFunc>()(v1 ^ c1.weight, v2 ^ c2.weight, v3 ^ c3.weight);
	}

template<typename Iterator>
	iterates_over_truth_table_t<Iterator>
	compute( node const& n, Iterator begin, Iterator end ) const
	{
		(void)end;

		assert( n != 0 && !is_ci( n ) );

		auto const& c1 = _storage->nodes[n].children[0];
		auto const& c2 = _storage->nodes[n].children[1];
		auto const& c3 = _storage->nodes[n].children[2];

		auto tt1 = *begin++;
		auto tt2 = *begin++;
		auto tt3 = *begin++;

		return compute_function<NodeFunc>()( c1.weight ? ~tt1 : tt1, c2.weight ? ~tt2 : tt2, c3.weight ? ~tt3 : tt3 );
	}
/*! \brief Re-compute the last block. */
template<typename Iterator>
	void compute( node const& n, kitty::partial_truth_table& result, Iterator begin, Iterator end ) const
	{
		static_assert( iterates_over_v<Iterator, kitty::partial_truth_table>, "begin and end have to iterate over partial_truth_tables" );

		(void)end;
		assert( n != 0 && !is_ci( n ) );

		auto const& c1 = _storage->nodes[n].children[0];
		auto const& c2 = _storage->nodes[n].children[1];
		auto const& c3 = _storage->nodes[n].children[2];

		auto tt1 = *begin++;
		auto tt2 = *begin++;
		auto tt3 = *begin++;

		assert( tt1.num_bits() > 0 && "truth tables must not be empty" );
		assert( tt1.num_bits() == tt2.num_bits() );
		assert( tt1.num_bits() == tt3.num_bits() );
		assert( tt1.num_bits() >= result.num_bits() );
		assert( result.num_blocks() == tt1.num_blocks() || ( result.num_blocks() == tt1.num_blocks() - 1 && result.num_bits() % 64 == 0 ) );

		result.resize( tt1.num_bits() );
		result._bits.back() = compute_function<NodeFunc>()( 
				c1.weight ? ~tt1._bits.back() : tt1._bits.back(),  
				c2.weight ? ~tt2._bits.back() : tt2._bits.back(),
				c3.weight ? ~tt3._bits.back() : tt3._bits.back());
		result.mask_bits();
	}

#pragma endregion

#pragma region Custom node values
				void clear_values() const
				{
					std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.data[0].h2 = 0; } );
				}

				auto value( node const& n ) const
				{
					return _storage->nodes[n].data[0].h2;
				}

				void set_value( node const& n, uint32_t v ) const
				{
					_storage->nodes[n].data[0].h2 = v;
				}

				auto incr_value( node const& n ) const
				{
					return _storage->nodes[n].data[0].h2++;
				}

				auto decr_value( node const& n ) const
				{
					return --_storage->nodes[n].data[0].h2;
				}
#pragma endregion

#pragma region Visited flags
				void clear_visited() const
				{
					std::for_each( _storage->nodes.begin(), _storage->nodes.end(), []( auto& n ) { n.data[1].h1 = 0; } );
				}

				auto visited( node const& n ) const
				{
					return _storage->nodes[n].data[1].h1;
				}

				void set_visited( node const& n, uint32_t v ) const
				{
					_storage->nodes[n].data[1].h1 = v;
				}

				uint32_t trav_id() const
				{
					return _storage->data.trav_id;
				}

				void incr_trav_id() const
				{
					++_storage->data.trav_id;
				}
#pragma endregion

#pragma region General methods
				auto& events() const
				{
					return *_events;
				}
#pragma endregion

			public:
				std::shared_ptr<tig_storage> _storage;
				std::shared_ptr<network_events<base_type>> _events;
		};

} // namespace mockturtle

namespace std
{

	template<>
		struct hash<mockturtle::tig_signal>
		{
			uint64_t operator()(mockturtle::tig_signal const& s ) const noexcept
			{
				uint64_t k = s.data;
				k ^= k >> 33;
				k *= 0xff51afd7ed558ccd;
				k ^= k >> 33;
				k *= 0xc4ceb9fe1a85ec53;
				k ^= k >> 33;
				return k;
			}
		}; /* hash */

	/*

	   Compile error: class template partial specialization contains a template parameter that cannot be deduced

	   using three_input_func =  bool (*)(bool, bool, bool);

	   template<three_input_function NodeFunc>
	   struct hash<mockturtle<NodeFunc>::tig_signal>
	   {
	   uint64_t operator()(mockturtle<NodeFunc>::tig_signal const& s ) const noexcept
	   {
	   uint64_t k = s.data;
	   k ^= k >> 33;
	   k *= 0xff51afd7ed558ccd;
	   k ^= k >> 33;
	   k *= 0xc4ceb9fe1a85ec53;
	   k ^= k >> 33;
	   return k;
	   }
	   }; 

*/


} // namespace std
