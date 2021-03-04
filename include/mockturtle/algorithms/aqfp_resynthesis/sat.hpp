#pragma once

#include <map>
#include <set>
#include <vector>

#include <bill/sat/interface/ghack.hpp>
#include <bill/sat/tseytin.hpp>
#include <fmt/format.h>

#include "./dag.hpp"

namespace mockturtle
{

using mig_chain = std::vector<std::pair<uint64_t, std::vector<uint32_t>>>;

uint64_t input_tt[] = {
    0x0000000000000000L,
    0xaaaaaaaaaaaaaaaaL,
    0xccccccccccccccccL,
    0xf0f0f0f0f0f0f0f0L,
    0xff00ff00ff00ff00L,
    0xffff0000ffff0000L,
    0xffffffff00000000L,
};

std::vector<uint64_t> primitives_3 = {
    0xe8UL,
    0xd4UL,
    0xb2UL,
    0x8eUL,
};

std::vector<uint64_t> primitives_5 = {
    0xfee8e880UL,
    0xfdd4d440UL,
    0xfbb2b220UL,
    0xf7717110UL,
    0xef8e8e08UL,
    0xdf4d4d04UL,
    0xbf2b2b02UL,
    0xe8fe80e8UL,
    0xd4fd40d4UL,
    0xb2fb20b2UL,
    0x8eef088eUL,
    0xe880fee8UL,
    0xd440fdd4UL,
    0xb220fbb2UL,
    0x8e08ef8eUL,
    0x80e8e8feUL,
};

std::map<uint32_t, std::vector<uint64_t>> primitives = {
    { 3u, primitives_3 },
    { 5u, primitives_5 },
};

std::vector<uint64_t> get_primitives( uint32_t fanin )
{
  assert( fanin == 3u || fanin == 5u );
  return primitives[fanin];
}

std::string as_string( bill::lbool_type val )
{
  switch ( val )
  {
  case bill::lbool_type::true_:
    return "1";
  case bill::lbool_type::false_:
    return "0";
  default:
    return "X";
  }
}

class sat_encoder
{

  void fix_normality()
  {
    output_inv = false;

    if ( ( tt & 1 ) != 0 )
    {
      tt = ( ~tt & ( ( 1UL << ( 1UL << num_in ) ) - 1 ) );
      output_inv = true;

      if ( verbose > 5 )
      {
        fmt::print( "function to synthesize is not normal. synthesizing the complement\n" );
      }
    }
  }

  void add_selection_variables()
  {
    for ( auto i = 0u; i <= num_in; i++ )
    {
      s_var.push_back( {} );
      for ( auto s = 0u; s < num_slots; s++ )
      {
        s_var[i].push_back(
            bill::lit_type( sol.add_variable(), bill::lit_type::polarities::positive ) );
      }
    }
  }

  void add_operator_variables()
  {
    for ( auto n = 0u; n < num_nodes; n++ )
    {
      f_var.push_back( {} );
      for ( auto j = 0u; j < ( 1 << net.node_num_fanin[n] ) - 1; j++ )
      {
        f_var[n].push_back(
            bill::lit_type( sol.add_variable(), bill::lit_type::polarities::positive ) );
      }
    }
  }

  void add_node_output_variables()
  {
    for ( auto n = 0u; n < num_all_nodes; n++ )
    {
      x_var.push_back( {} );
      for ( auto r = 0u; r < num_rows; r++ )
      {
        x_var[n].push_back(
            bill::lit_type( sol.add_variable(), bill::lit_type::polarities::positive ) );
      }
    }
  }

  void add_selection_clauses()
  {
    // no input is connected to more than one slot
    if ( !net.is_partial_dag )
    {
      for ( auto i = 0u; i <= num_in; i++ )
      {
        for ( auto s1 = 0u; s1 < num_slots; s1++ )
        {
          for ( auto s2 = s1 + 1; s2 < num_slots; s2++ )
          {
            clauses.push_back(
                ~add_tseytin_and( sol, s_var[i][s1], s_var[i][s2] ) );
          }
        }
      }
    }

    if ( !net.is_partial_dag )
    {
      if ( net.zero_input == 0 )
      { // zero input not connected anywhere
        std::vector<bill::lit_type> temp;
        for ( auto s = 0u; s < num_slots; s++ )
        {
          temp.push_back( ~s_var[0][s] );
        }
        clauses.push_back( add_tseytin_and( sol, temp ) );
      }
      else
      { // zero input connected to net.zero_input
        clauses.push_back( s_var[0][net.zero_input - net.num_gates()] );
      }
    }

    // each slot connected to one input
    for ( auto s = 0u; s < num_slots; s++ )
    {
      std::vector<bill::lit_type> temp;
      for ( auto i = 0u; i <= num_in; i++ )
      {
        temp.push_back( s_var[i][s] );
      }
      clauses.push_back( add_tseytin_or( sol, temp ) );
    }

    // no slot is connected to more than one input
    for ( auto i1 = 0u; i1 <= num_in; i1++ )
    {
      for ( auto i2 = i1 + 1u; i2 <= num_in; i2++ )
      {
        for ( auto s = 0u; s < num_slots; s++ )
        {
          clauses.push_back(
              ~add_tseytin_and( sol, s_var[i1][s], s_var[i2][s] ) );
        }
      }
    }
  }

  void add_input_truthtable_clauses()
  {
    // clauses for truthtables of input slots
    for ( auto i = 0u; i <= num_in; i++ )
    {
      for ( auto s = 0u; s < num_slots; s++ )
      {
        for ( auto r = 0u; r < num_rows; r++ )
        {
          // s[i][s] => ( x[s,r] = input_tt[tt_ind][r+1] )
          if ( ( input_tt[i] >> ( r + 1 ) ) & 1 )
          {
            clauses.push_back(
                add_tseytin_or( sol, ~s_var[i][s],
                                x_var[s + num_nodes][r] ) );
          }
          else
          {
            clauses.push_back(
                add_tseytin_or( sol, ~s_var[i][s],
                                ~x_var[s + num_nodes][r] ) );
          }
        }
      }
    }
  }

  void add_node_operator_clauses()
  {
    for ( auto n = 0u; n < num_nodes; n++ )
    {
      std::vector<bill::lit_type> temp2;
      for ( auto p : get_primitives( net.node_num_fanin[n] ) )
      {
        std::vector<bill::lit_type> temp;
        for ( auto j = 0u; j < ( 1 << net.node_num_fanin[n] ) - 1; j++ )
        {
          if ( ( p >> ( j + 1 ) ) & 1 )
          {
            temp.push_back( f_var[n][j] );
          }
          else
          {
            temp.push_back( ~f_var[n][j] );
          }
        }
        temp2.push_back( add_tseytin_and( sol, temp ) );
      }

      clauses.push_back( add_tseytin_or( sol, temp2 ) );
    }
  }

  void add_node_input_output_operator_consistency_clauses()
  {
    for ( auto n = 0u; n < num_nodes; n++ )
    {
      assert( net.node_num_fanin[n] == 3u || net.node_num_fanin[n] == 5u );

      if ( net.node_num_fanin[n] == 3u )
      {
        for ( auto r = 0u; r < num_rows; r++ )
        {
          const auto t0 = x_var[net.nodes[n][0]][r];
          const auto t1 = x_var[net.nodes[n][1]][r];
          const auto t2 = x_var[net.nodes[n][2]][r];
          const auto t3 = x_var[n][r];

          std::vector<bill::lit_type> temp = { t0, t1, t2, ~t3 };

          clauses.push_back( add_tseytin_or( sol, temp ) );

          for ( auto d = 0u; d < 2u; d++ )
          {
            for ( auto c = 0u; c < 2u; c++ )
            {
              for ( auto b = 0u; b < 2u; b++ )
              {
                for ( auto a = 0u; a < 2u; a++ )
                {
                  auto f_ind = a + ( 2 * b ) + ( 4 * c );

                  if ( f_ind == 0 )
                  {
                    continue;
                  }
                  f_ind--;

                  // n = (n0, n1, n2)
                  // (x_{n0,r} = a) and (x_{n1,r} = b) and
                  // (x_{n2,r} = c) and (x_{n,r} = d) => (f_{c,b,a} = d)

                  const auto t4 = f_var[n][f_ind];

                  std::vector<bill::lit_type> temp2 = { ( a == 0 ) ? t0 : ~t0,
                                                        ( b == 0 ) ? t1 : ~t1,
                                                        ( c == 0 ) ? t2 : ~t2,
                                                        ( d == 1 ) ? t3 : ~t3,
                                                        ( d == 0 ) ? t4 : ~t4 };

                  clauses.push_back(
                      add_tseytin_or( sol, temp2 ) );
                }
              }
            }
          }
        }
      }
      else
      {
        for ( auto r = 0u; r < num_rows; r++ )
        {
          const auto t0 = x_var[net.nodes[n][0]][r];
          const auto t1 = x_var[net.nodes[n][1]][r];
          const auto t2 = x_var[net.nodes[n][2]][r];
          const auto t3 = x_var[net.nodes[n][3]][r];
          const auto t4 = x_var[net.nodes[n][4]][r];
          const auto t5 = x_var[n][r];

          std::vector<bill::lit_type> temp = { t0, t1, t2, t3, t4, ~t5 };

          clauses.push_back( add_tseytin_or( sol, temp ) );

          for ( auto f = 0u; f < 2u; f++ )
          {
            for ( auto e = 0u; e < 2u; e++ )
            {
              for ( auto d = 0u; d < 2u; d++ )
              {
                for ( auto c = 0u; c < 2u; c++ )
                {
                  for ( auto b = 0u; b < 2u; b++ )
                  {
                    for ( auto a = 0u; a < 2u; a++ )
                    {
                      auto f_ind = a + ( 2u * b ) + ( 4u * c ) + ( 8u * d ) + ( 16u * e );

                      if ( f_ind == 0 )
                      {
                        continue;
                      }
                      f_ind--;

                      const auto t6 = f_var[n][f_ind];

                      std::vector<bill::lit_type> temp2 = { ( a == 0 ) ? t0 : ~t0,
                                                            ( b == 0 ) ? t1 : ~t1,
                                                            ( c == 0 ) ? t2 : ~t2,
                                                            ( d == 0 ) ? t3 : ~t3,
                                                            ( e == 0 ) ? t4 : ~t4,
                                                            ( f == 1 ) ? t5 : ~t5,
                                                            ( f == 0 ) ? t6 : ~t6 };

                      clauses.push_back(
                          add_tseytin_or( sol, temp2 ) );
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  void add_network_output_clauses()
  {
    for ( auto r = 0u; r < num_rows; r++ )
    {
      // (i + 1)-th row of top nodes truthtable must equal that of tt
      if ( ( tt >> ( r + 1 ) ) & 1 )
      {
        clauses.push_back( x_var[0][r] );
      }
      else
      {
        clauses.push_back( ~x_var[0][r] );
      }
    }
  }

public:
  sat_encoder( const aqfp_dag<>& _net, uint32_t _num_in, uint64_t _tt, uint32_t _verbose = 0 )
  {
    net = _net;
    num_in = _num_in;
    tt = _tt;
    verbose = _verbose;

    num_slots = net.input_slots.size();
    num_nodes = net.num_gates();
    num_all_nodes = net.nodes.size();
    num_rows = ( 1 << num_in ) - 1;

    fix_normality();
    if ( verbose > 5 )
    {
      fmt::print( "synthesizing function: {:08x}\n", tt );
    }

    /* NOTE: do not change the order of following function calls */

    add_selection_variables();
    if ( verbose > 5 )
    {
      fmt::print( "added selection variables\n" );
    }

    add_operator_variables();
    if ( verbose > 5 )
    {
      fmt::print( "added operator variables\n" );
    }

    add_node_output_variables();
    if ( verbose > 5 )
    {
      fmt::print( "added node output variables\n" );
    }

    add_selection_clauses();
    if ( verbose > 5 )
    {
      fmt::print( "added selection clauses\n" );
    }

    add_input_truthtable_clauses();

    if ( verbose > 5 )
    {
      fmt::print( "added input truthtable clauses\n" );
    }

    add_node_operator_clauses();
    if ( verbose > 5 )
    {
      fmt::print( "added node operator clauses\n" );
    }

    add_node_input_output_operator_consistency_clauses();
    if ( verbose > 5 )
    {
      fmt::print( "added node input, output, and operator consistency clauses\n" );
    }

    add_network_output_clauses();
    if ( verbose > 5 )
    {
      fmt::print( "added network output clauses\n" );
    }

    const auto root = add_tseytin_and( sol, clauses );

    sol.add_clause( root );
  }

  bool solve()
  {
    if ( sol.solve() == bill::result::states::unsatisfiable )
    {
      return false;
    }
    else
    {
      return true;
    }
  }

  std::string chain_as_string()
  {
    std::map<uint32_t, uint32_t> input_map;
    std::map<uint32_t, uint64_t> node_func_map;

    auto model = sol.get_result().model();
    uint32_t var_ind = 0;

    for ( auto i = 0u; i <= num_in; i++ )
    {
      for ( auto s = 0u; s < num_slots; s++ )
      {
        if ( model[var_ind++] == bill::lbool_type::true_ )
        {
          input_map[num_nodes + s] = i;
        }
      }
    }

    for ( auto n = 0u; n < num_nodes; n++ )
    {
      uint64_t val = 0;
      for ( auto j = 0u; j < ( 1 << net.node_num_fanin[n] ) - 1; j++ )
      {
        val |= ( 1u << ( j + 1 ) ) * ( ( model[var_ind++] == bill::lbool_type::true_ ) ? 1 : 0 );
      }
      node_func_map[n] = val;
    }

    mig_chain chain( num_nodes );
    for ( auto n = 0u; n < num_nodes; n++ )
    {
      chain[num_nodes - 1 - n].first = node_func_map[n];

      for ( auto i = 0u; i < net.node_num_fanin[n]; i++ )
      {
        if ( (uint32_t)net.nodes[n][i] >= num_nodes )
        {
          chain[num_nodes - 1 - n].second.push_back(
              input_map[net.nodes[n][i]] );
        }
        else
        {
          chain[num_nodes - 1 - n].second.push_back(
              num_nodes - 1 - net.nodes[n][i] + num_in + 1 );
        }
      }
    }
    std::string res = "";
    res += fmt::format( "{} {}\n", num_in, num_nodes );
    for ( auto&& c : chain )
    {
      res += fmt::format( "{:8x} {}\n", c.first, fmt::join( c.second, " " ) );
    }
    return res;
  }

  void print_solution()
  {
    auto model = sol.get_result().model();
    uint32_t var_ind = 0;
    for ( auto i = 0u; i <= num_in; i++ )
    {
      fmt::print( "i = {} x_is =", i );
      for ( auto s = 0u; s < num_slots; s++ )
      {
        fmt::print( " {}", as_string( model[var_ind++] ) );
      }
      fmt::print( "\n" );
    }

    for ( auto n = 0u; n < num_nodes; n++ )
    {
      std::string tmp = "0";
      uint32_t val = 0;
      for ( auto j = 0u; j < 7u; j++ )
      {
        val |= ( 1 << ( j + 1 ) ) * ( ( model[var_ind] == bill::lbool_type::true_ ) ? 1 : 0 );
        tmp = as_string( model[var_ind++] ) + tmp;
      }
      fmt::print( "n = {} f_n = {} ({})\n", n, tmp, val );
    }

    for ( auto n = 0u; n < num_all_nodes; n++ )
    {
      std::string tmp = "0";
      for ( auto r = 0u; r < num_rows; r++ )
      {
        tmp = as_string( model[var_ind++] ) + tmp;
      }
      fmt::print( "tt n = {} tt = {}\n", n, tmp );
    }
  }

  uint64_t get_tt()
  {
    return tt;
  }

private:
  aqfp_dag<> net;

  bool output_inv;
  uint64_t tt;

  uint32_t verbose;

  uint32_t num_in;
  uint32_t num_slots;
  uint32_t num_nodes;
  uint32_t num_all_nodes;
  uint32_t num_rows;

  std::vector<bill::lit_type> clauses;

  std::vector<std::vector<bill::lit_type>> s_var; // which variable or constant is connected to which input
  std::vector<std::vector<bill::lit_type>> f_var; // which gate computes which primitive
  std::vector<std::vector<bill::lit_type>> x_var; // output truthtables

  bill::solver<bill::solvers::ghack> sol;
};

} // namespace mockturtle
