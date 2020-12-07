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
  \file exact.hpp
  \brief Replace with exact synthesis result

  \author Mathias Soeken
*/

#pragma once

#include <iostream>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/hash.hpp>
#include <kitty/print.hpp>
#include <kitty/traits.hpp>

#include "../../networks/tig.hpp"
#include "../../networks/tig_and3.hpp"
#include "../../networks/tig_andxor.hpp"
#include "../../networks/tig_dot.hpp"
#include "../../networks/tig_gamble.hpp"
#include "../../networks/tig_maj.hpp"
#include "../../networks/tig_mux.hpp"
#include "../../networks/tig_onehot.hpp"
#include "../../networks/tig_orand.hpp"
#include "../../networks/tig_xorand.hpp"

#include "../../utils/include/percy.hpp"

namespace mockturtle
{

template<three_input_function NodeFunc, typename Ntk>
class exact_tig_resynthesis
{
public:
  explicit exact_tig_resynthesis( exact_tig_resynthesis_params const& ps = {} ) : ps( ps )
  {
	  prim_map = std::map<uint64_t, std::tuple<int, int, int, bool, kitty::dynamic_truth_table>>();

	  kitty::dynamic_truth_table cnst0{3};
	  kitty::dynamic_truth_table a{3};
	  kitty::dynamic_truth_table b{3};
	  kitty::dynamic_truth_table c{3};
	  kitty::create_nth_var( a, 0 );
	  kitty::create_nth_var( b, 1 );
	  kitty::create_nth_var( c, 2 );

	  std::array<const kitty::dynamic_truth_table, 8> choices = {a, ~a, b, ~b, c, ~c, cnst0, ~cnst0};

	  for (auto ind0 = 0; ind0 < 8; ind0++)
	  {
		  for (auto ind1 = 0; ind1 < 8; ind1++)
		  {
			  for (auto ind2 = 0; ind2 < 8; ind2++)
			  { 
				  kitty::dynamic_truth_table prim = compute_function<NodeFunc>()(choices[ind0], choices[ind1], choices[ind2]); 
				  bool negated = false;
				  if (!kitty::is_normal(prim)) {
					  prim = ~prim;
					  negated = true;
				  }
				  prim_map[prim._bits[0]] = {ind0, ind1, ind2, negated, prim};
			  }
		  }
	  }	
  }

  template<typename LeavesIterator, typename TT, typename Fn>
  void operator()( Ntk& ntk, TT const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn ) const
  {
	  static_assert( kitty::is_complete_truth_table<TT>::value, "Truth table must be complete" );

	  using signal = mockturtle::signal<Ntk>;

	  auto const tt = function.num_vars() < 4u ? kitty::extend_to( function, 4u ) : function;

	  percy::chain chain;
	  bool found_in_cache = false;

	  const auto config = kitty::exact_npn_canonization( tt );
	  auto npn = std::get<0>(config);

	  bool const normal = kitty::is_normal( npn );

	  if ( ps.cache ) {
		  if (const auto it = ps.cache->find( normal ? npn : ~npn ); it != ps.cache->end() ) {
			  chain = it->second;
			  found_in_cache = true;
		  } 
	  }

	  if (!found_in_cache) {
		  percy::spec spec;
		  spec.verbosity = 0;
		  spec.fanin = 3;
		  spec[0] = normal ? npn : ~npn;

		  for (auto it = prim_map.begin(); it != prim_map.end(); it++) {
			  spec.add_primitive( std::get<4>(it->second) );
		  }	

		  const auto result = percy::synthesize( spec, chain, ps.solver_type, ps.encoder_type, ps.synthesis_method );
		  assert( result == percy::success );

		  auto const sim = chain.simulate();
		  assert( chain.simulate()[0] == spec[0] );
	  }

	  std::vector<signal> pis( 4, ntk.get_constant( false ) );
	  std::copy( begin, end, pis.begin() );

	  std::vector<signal> signals( 4, ntk.get_constant( false ) );
	  auto perm = std::get<2>( config );
	  for ( auto i = 0u; i < 4; ++i )
	  {
		  signals[i] = pis[perm[i]];
	  }

	  const auto& phase = std::get<1>( config );
	  for ( auto i = 0u; i < 4; ++i )
	  {
		  if ( ( phase >> perm[i] ) & 1 )
		  {
			  signals[i] = !signals[i];
		  }
	  }

	  for ( auto i = 0; i < chain.get_nr_steps(); ++i )
	  {
		  auto const c1 = signals[chain.get_step( i )[0]];
		  auto const c2 = signals[chain.get_step( i )[1]];
		  auto const c3 = signals[chain.get_step( i )[2]];
		  
		  signal sigs[8] = { c1, !c1, c2, !c2, c3, !c3, ntk.get_constant( false ), !ntk.get_constant( false ) };

		  uint64_t op = chain.get_operator( i )._bits[0];
		  if (prim_map.count(op)) {
		  	  std::tuple<int, int, int, bool, kitty::dynamic_truth_table> tup = prim_map.find(op)->second;
			  auto i1 = std::get<0>(tup);
			  auto i2 = std::get<1>(tup);
			  auto i3 = std::get<2>(tup);
			  auto inv = std::get<3>(tup);

			  auto sig = ntk.create_network_node( sigs[i1], sigs[i2], sigs[i3] );

			  signals.emplace_back( inv ? !sig : sig );
		  } else {
			  fmt::print("need to construct {} but not found!\n", op);
			  for (auto it = prim_map.begin(); it != prim_map.end(); it++) {
				  fmt::print("pattern {}\n", it->first);
			  }
			  assert( false );
		  }
	  }

	  assert( chain.get_outputs().size() > 0u );
	  uint32_t const output_index = ( chain.get_outputs()[0u] >> 1u );
	  auto const output_signal = output_index == 0u ? ntk.get_constant( false ) : signals[output_index - 1];

	  fn( chain.is_output_inverted( 0 ) ^ normal  ^ ( ( phase >> tt.num_vars() ) & 1 ) ? output_signal : !output_signal ); 
  }

  std::map<uint64_t, std::tuple<int, int, int, bool, kitty::dynamic_truth_table>> prim_map;

protected:
  exact_tig_resynthesis_params const& ps;
}; /* exact_dig_resynthesis */

} /* namespace mockturtle */

