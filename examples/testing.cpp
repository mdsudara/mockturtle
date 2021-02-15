#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>

#include <fmt/format.h>

#include <lorina/lorina.hpp>

#include <mockturtle/algorithms/exact_syn/aqfp_db.hpp>
#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <mockturtle/algorithms/exact_syn/dag_cost.hpp>
#include <mockturtle/algorithms/exact_syn/gen_dag.hpp>
#include <mockturtle/algorithms/exact_syn/generate_db.hpp>
#include <mockturtle/algorithms/exact_syn/sat.hpp>
#include <mockturtle/mockturtle.hpp>

#include "utils.hpp"

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::set<T>& s );

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::unordered_set<T>& s );

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::multiset<T>& s );

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::vector<T>& s );

template<typename T, typename S>
std::ostream& operator<<( std::ostream& os, const std::map<T, S>& s );

template<typename T, typename S>
std::ostream& operator<<( std::ostream& os, const std::unordered_map<T, S>& s );

template<typename T>
std::ostream& operator<<( std::ostream& os, const mockturtle::aqfp_logical_network_t<T>& s );

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::set<T>& s )
{
  os << "{";
  bool first = true;
  for ( const T& t : s )
  {
    if ( not first )
      os << ",";
    os << " " << t;
    first = false;
  }
  os << " }";

  return os;
}

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::unordered_set<T>& s )
{
  os << "{";
  bool first = true;
  for ( const T& t : s )
  {
    if ( not first )
      os << ",";
    os << " " << t;
    first = false;
  }
  os << " }";

  return os;
}

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::multiset<T>& s )
{
  os << "{";
  bool first = true;
  for ( const T& t : s )
  {
    if ( not first )
      os << ",";
    os << " " << t;
    first = false;
  }
  os << " }";

  return os;
}

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::vector<T>& s )
{
  os << "{";
  bool first = true;
  for ( const T& t : s )
  {
    if ( not first )
      os << ",";
    os << " " << t;
    first = false;
  }
  os << " }";

  return os;
}

template<typename T, typename S>
std::ostream& operator<<( std::ostream& os, const std::map<T, S>& s )
{
  os << "{";
  bool first = true;
  for ( auto it = s.begin(); it != s.end(); it++ )
  {
    if ( not first )
      os << ",";
    os << " "
       << "{ " << it->first << ", " << it->second << " }";
    first = false;
  }
  os << " }";

  return os;
}

template<typename T, typename S>
std::ostream& operator<<( std::ostream& os, const std::unordered_map<T, S>& s )
{
  os << "{";
  bool first = true;
  for ( auto it = s.begin(); it != s.end(); it++ )
  {
    if ( not first )
      os << ",";
    os << " "
       << "{ " << it->first << ", " << it->second << " }";
    first = false;
  }
  os << " }";

  return os;
}

template<typename T>
std::ostream& operator<<( std::ostream& os, const mockturtle::aqfp_logical_network_t<T>& s )
{
  os << "{ ";
  os << s.node_num_fanin << ", ";
  os << s.nodes << ", ";
  os << s.input_slots << ", ";
  os << s.zero_input << ", ";
  os << s.pdag_id << "u, ";
  os << s.is_partial_dag << ", ";
  os << s.num_levels << "u, ";
  os << s.computed_cost << ", ";
  os << s.num_gates_of_fanin << ", ";
  os << s.last_layer_leaves << ", ";
  os << s.other_leaves << "}";

  return os;
}

void test_partition_generation()
{
  mockturtle::detail::partition_generator<int> partition_gen;

  std::cout << partition_gen( { 0, 1, 2, 3 }, { 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 2, 3 }, { 1u, 1u, 1u, 1u }, 0u, 0u ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 0u ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 2u ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 0u, 2u ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 1, 3 }, { 1u, 1u, 1u, 1u }, 2u, 2u ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 1, 3 }, { 1u, 2u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_gen( { 0, 1, 1, 3 }, { 1u, 2u, 1u, 1u }, 0u, 2u ) << std::endl
            << std::endl;
}

void test_partition_extension()
{
  mockturtle::detail::partition_extender<int> partition_ext;

  std::cout << partition_ext( { 0, 1, 2, 3, 4 }, {}, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( {}, { { 1, 2 }, { 3, 4 } }, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( { 0, 1, 2, 3, 4 }, { {} }, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( { 4 }, { { 1 }, { 2 }, { 3 }, { 4 } }, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( { 4 }, { { 1 }, { 2 }, { 3 }, { 4 } }, { 1u, 1u, 1u, 1u, 2u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( { 3, 4 }, { { 0 }, { 1 }, { 2 } }, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( { 0, 1, 2 }, { { 3 }, { 4 } }, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
  std::cout << partition_ext( { 3, 4, 4, 4 }, { { 0 }, { 1 }, { 2 } }, { 1u, 1u, 1u, 1u, 1u } ) << std::endl
            << std::endl;
}

void test_sublist_generation()
{
  mockturtle::detail::sublist_generator<int> sublist_gen;
  std::cout << sublist_gen( { 0 } ) << std::endl
            << std::endl;
  std::cout << sublist_gen( { 0, 0, 0, 0 } ) << std::endl
            << std::endl;
  std::cout << sublist_gen( { 0, 1, 2, 3 } ) << std::endl
            << std::endl;
  std::cout << sublist_gen( { 0, 1, 1, 3 } ) << std::endl
            << std::endl;
  std::cout << sublist_gen( { 0, 2, 2, 2 } ) << std::endl
            << std::endl;
  std::cout << sublist_gen( { 0, 0, 2, 2, 2 } ) << std::endl
            << std::endl;
}

void test_layer_extension()
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  std::vector<uint32_t> allowed_num_fanins = { 3u, 5u };
  std::map<uint32_t, uint32_t> max_gates_of_fanin = { { 3u, 3u }, { 5u, 1u } };

  auto params = mockturtle::dag_generator_params();

  params.max_gates = 3u;         // allow at most 4 gates in total
  params.max_num_fanout = 1000u; // limit the maximum fanout of a gate
  params.max_width = 1000u;      // maximum number of gates at any level
  params.max_num_in = 3u;        // maximum number of inputs slots (need extra one for the constant)
  params.max_levels = 3u;        // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  auto dummy_cost_computer = []( const auto& net ) {(void)net; return 0.0; };

  mockturtle::dag_generator<int, decltype( dummy_cost_computer )> gen( params, dummy_cost_computer );
  Ntk pdag1 = {
      { 3, 3 },
      { { 1 }, {} },
      {},
      0,
      0,
      true,
      2u,
      -1.0,
      { { 3u, 2u } },
      { 1, 1, 1 },
      { 0, 0 } };

  std::cout << "Layer extension 1\n\n";
  for ( auto x : gen.get_layer_extension( pdag1 ) )
  {
    std::cout << x << "\n";
  }

  Ntk pdag2 = {
      { 5 },
      { {} },
      {},
      0,
      0,
      true,
      1u,
      -1.0,
      { { 5u, 1u } },
      { 0, 0, 0, 0, 0 },
      {} };

  std::cout << "\n\nLayer extension 2\n\n";
  for ( auto x : gen.get_layer_extension( pdag2 ) )
  {
    std::cout << x << "\n";
  }

  std::cout << "\n";

  std::cout << "\n\nDAG generation 1\n\n";
  for ( auto x : gen.get_dags_from_partial_dag( pdag1 ) )
  {
    std::cout << x << "\n";
  }

  std::cout << "\n\nDAG generation 2\n\n";
  for ( auto x : gen.get_dags_from_partial_dag( pdag2 ) )
  {
    std::cout << x << "\n";
  }

  std::cout << "\n";
}

void test_dag_simulation()
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;
  Ntk dag1 = { { 3, 3, 0, 0, 0 }, { { 1, 2, 3 }, { 2, 3, 4 }, {}, {}, {} }, { 2, 3, 4 }, 0, 0u, 0, 2u, -1, { { 0, 3 }, { 3, 2 } }, {}, {} };
  Ntk dag2 = { { 3, 3, 0, 0, 0 }, { { 1, 2, 3 }, { 2, 3, 4 }, {}, {}, {} }, { 2, 3, 4 }, 3, 0u, 0, 2u, -1, { { 0, 3 }, { 3, 2 } }, {}, {} };

  mockturtle::dag_simulator<uint64_t> sim( { 0x00ul, 0xaaul, 0xccul, 0xf0ul } );

  std::cout << sim.all_functions_from_dag( dag1 ) << std::endl;
  std::cout << sim.all_functions_from_dag( dag2 ) << std::endl;
}

void test_timing()
{
  using Ntk = mockturtle::aqfp_logical_network_t<int>;

  std::vector<uint32_t> allowed_num_fanins = { 3u };
  std::map<uint32_t, uint32_t> max_gates_of_fanin = { { 3u, 7u } };

  auto params = mockturtle::dag_generator_params();

  params.max_gates = 6u;  // allow at most 4 gates in total
  params.max_num_in = 4u; // maximum number of inputs slots (need extra one for the constant)
  params.max_levels = 3u; // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  auto trivial_cost = []( auto& net ) {(void) net; return 0.0; };
  mockturtle::dag_generator<int, decltype( trivial_cost )> gen( params, trivial_cost );

  mockturtle::aqfp_cost_computer<Ntk> aqfp_cc( { { 3u, 3.0 } }, { { 3u, 3.0 } }, 1.0, 3u );

  std::vector<uint64_t> input_tt = {
      0x0000UL,
      0xaaaaUL,
      0xccccUL,
      0xf0f0UL,
      0xff00UL,
  };
  mockturtle::dag_simulator<uint64_t> dag_sim( input_tt );

  std::vector<Ntk> all_dags;

  auto t0 = std::chrono::high_resolution_clock::now();

  auto c0 = 0u;
  while ( true )
  {
    auto x = gen.next_dag();
    if ( x == std::nullopt )
      break;

    c0++;
    if ( x.value().num_gates() == 6u )
    {
      all_dags.push_back( x.value() );
    }
  }

  auto t1 = std::chrono::high_resolution_clock::now();

  auto c1 = 0u;
  for ( auto& dag : all_dags )
  {
    aqfp_cc( dag );
    if ( ( ++c1 ) >= 1000u )
      break;
  }

  auto t2 = std::chrono::high_resolution_clock::now();

  auto c2 = 0u;
  for ( auto& dag : all_dags )
  {
    dag_sim.all_functions_from_dag( dag );
    if ( ( ++c2 ) >= 1000u )
      break;
  }

  auto t3 = std::chrono::high_resolution_clock::now();

  std::vector<uint64_t> tts = { 0x03fc, 0x19e6, 0x06b4, 0xfc26, 0x0368 };
  auto c3 = 0u;
  for ( auto& dag : all_dags )
  {
    for ( auto tt : tts )
    {
      mockturtle::sat_encoder se( dag, 4u, tt, 0u );
      se.solve();
    }
    if ( ( ++c3 ) >= 20u )
      break;
  }

  auto t4 = std::chrono::high_resolution_clock::now();

  auto d1 = std::chrono::duration_cast<std::chrono::microseconds>( t1 - t0 ).count() / 1e6;
  auto d2 = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() / 1e6;
  auto d3 = std::chrono::duration_cast<std::chrono::microseconds>( t3 - t2 ).count() / 1e6;
  auto d4 = std::chrono::duration_cast<std::chrono::microseconds>( t4 - t3 ).count() / 1e6;

  fmt::print( "Generating {:5d} DAGs took                   {:5.4f} seconds\n", c0, d1 );
  fmt::print( "Computing costs for {:4d} DAGs took           {:5.4f} seconds\n", c1, d2 );
  fmt::print( "Simulating {:4d} DAGs took                    {:5.4f} seconds\n", c2, d3 );
  fmt::print( "Synthesizing {:1d} functions on {:2d} DAGs took     {:5.4f} seconds\n", tts.size(), c3, d4 );
}

void assemble_db()
{
  std::string dags_file_prefix = "/Volumes/Dewmini_5/AQFP/Dags/dags_";
  std::string cost_file_prefix = "/Volumes/Dewmini_5/AQFP/Costs-G6-S42-B2/cost_dags_";

  std::vector<std::string> suff = { "aa", "ab", "ac", "ad", "ae", "af", "ag", "ah",
                                    "ai", "aj", "ak", "al", "am", "an", "ao", "ap",
                                    "aq", "ar", "as", "at", "au", "av", "aw", "ax",
                                    "ay", "az", "ba", "bb", "bc", "bd", "be", "bf",
                                    "bg", "bh", "bi", "bj", "bk", "bl", "bm", "bn",
                                    "bo", "bp", "bq", "br", "bs", "bt", "bu", "bv" };

  //  std::vector<std::string> suff = { "aa", "ab", "ac", "ad", "ae"};

  mockturtle::aqfp_logical_network_db<mockturtle::aqfp_logical_network_t<int>, 4u> db( 1u );

  for ( std::string suf : suff )
  {
    std::ifstream d( dags_file_prefix + suf );
    std::ifstream c( cost_file_prefix + suf );

    fmt::print( "{} {}\n", d.is_open(), c.is_open() );

    std::string dag;
    auto dag_count = 0u;

    while ( std::getline( d, dag ) )
    {
      dag_count++;
      auto ok = dag_count < 100 || dag_count % 10023 == 0u;

      if ( ok )
      {
        std::string token;
        c >> token;
        std::unordered_map<uint64_t, double> configs;
        while ( true )
        {
          double cst;
          c >> token;
          if ( token == "end" )
          {
            break;
          }
          c >> cst;
          configs[std::stoul( token, 0, 16 )] = cst;
        }

        mockturtle::aqfp_logical_network_t<int> ntk;
        ntk.decode_dag( dag );
        if ( ntk.input_slots.size() < 5u || ( ntk.input_slots.size() == 5u && ntk.zero_input != 0 ) )
        {
          db.update( ntk, configs );
        }
      }
      else
      {
        std::string token;
        c >> token;
        while ( true )
        {
          c >> token;
          if ( token == "end" )
          {
            break;
          }
          c >> token;
        }
      }

      if ( dag_count % 100000 == 0u )
      {
        std::cerr << fmt::format( "count = {}\n", dag_count );
      }

      if ( dag_count == 10000000 )
      {
        break;
      }
    }
    std::string outf = "db_" + suf + ".txt";
    std::ofstream dbout( outf );
    db.save_db_to_file( dbout );
    dbout.close();

    d.close();
    c.close();
  }
}

void reverse_permutations() {
  std::vector<std::string> suff = {
      "aa", "ab", "ac", "ad", "ae", "af", "ag", "ah",
      "ai", "aj", "ak", "al", "am", "an", "ao", "ap",
      "aq", "ar", "as", "at", "au", "av", "aw", "ax",
      "ay", "az", "ba", "bb", "bc", "bd", "be", "bf",
      "bg", "bh", "bi", "bj", "bk", "bl", "bm", "bn",
      "bo", "bp", "bq", "br", "bs", "bt", "bu", "bv",
  };


  for ( auto suf : suff )
  {
    mockturtle::aqfp_logical_network_db<mockturtle::aqfp_logical_network_t<int>, 4u> db( 2.0, 1u );

    std::string inf = "/Users/MDS/EPFL/repos/FromServer/DB-complete/db_dags_" + suf;
    fmt::print( "reading db {}\n", inf );

    std::ifstream dbchunk( inf );
    assert( dbchunk.is_open() );
    db.load_db_from_file( dbchunk );
    dbchunk.close();

    std::ofstream final_db( "db_dags_" + suf + "_inputs_reversed.txt" );
    assert(final_db.is_open());
    db.save_db_to_file( final_db );
    final_db.close();;
  }
}

void combine_db()
{
  std::vector<std::string> suff = {
       "aa", "ab", "ac", "ad", "ae", "af", "ag", "ah",
       "ai", "aj", "ak", "al", "am", "an", "ao", "ap",
       "aq", "ar", "as", "at", "au", "av", "aw", "ax",
       "ay", "az", "ba", "bb", "bc", "bd", "be", "bf",
       "bg", "bh", "bi", "bj", "bk", "bl", "bm", "bn",
       "bo", "bp", "bq", "br", "bs", "bt", "bu", "bv",
  };

  mockturtle::aqfp_logical_network_db<mockturtle::aqfp_logical_network_t<int>, 4u> db( 2.0, 1u );

  for ( auto suf : suff )
  {
    std::string inf = "/Users/dewmini/mockturtle/DB/db_dags_" + suf;
    // std::string inf = "db_dags_" + suf + "_inputs_reversed.txt";
    fmt::print( "reading db {}\n", inf );
    std::ifstream dbchunk( inf );
    assert( dbchunk.is_open() );
    db.load_db_from_file( dbchunk );
    dbchunk.close();
  }

  fmt::print( "removing redundancies\n" );
  db.remove_redundant();

  std::ofstream final_db( "aqfp_db_2020_feb_05.txt" );
  db.save_db_to_file( final_db );
  final_db.close();
}

uint32_t num_splitter_levels( uint32_t fanouts, uint32_t splitter_cap )
{
  return std::ceil( std::log( fanouts ) / std::log( splitter_cap ) );
}

uint32_t num_splitters( uint32_t fanouts, uint32_t splitter_cap )
{
  //  1 + ( splitter_cap - 1 ) * k >= fanouts;
  return ( fanouts - 1 + splitter_cap - 2 ) / ( splitter_cap - 1 );
}

// void test_aqfp_syn()
// {
//   mockturtle::mig_network ntk;
//   lorina::read_verilog( "/Users/MDS/EPFL/repos/Benchmarks/bar.v", mockturtle::verilog_reader( ntk ) );
//   fmt::print( "pi: {} po: {} size: {}\n", ntk.num_pis(), ntk.num_pos(), ntk.num_gates() );
//   aqfp_view aqfp( ntk );
//   fmt::print( "" );

//   mockturtle::klut_network klut = lut_map_old( ntk, "ctrl", "rewrite", 4 );
//   mockturtle::depth_view lutd{ klut };
//   fmt::print( "depth {} \n", lutd.depth() );

//   mockturtle::aqfp_logical_network_db<mockturtle::aqfp_logical_network_t<int>, 4u> db( 1u );
//   std::ifstream final_db( "aqfp_db.txt" );
//   assert( final_db.is_open() );
//   db.load_db_from_file( final_db );
//   final_db.close();

//   // std::ofstream final_db2( "aqfp_db_compressed2.txt");
//   // db.save_db_to_file(final_db2);
//   // final_db2.close();

//   double total_cost = 0;

//   std::unordered_map<uint64_t, uint32_t> level;
//   std::unordered_map<uint64_t, uint32_t> fanout;

//   klut.foreach_pi( [&level]( auto n ) { level[n] = 0; } );

//   std::vector<uint64_t> nodes;
//   klut.foreach_node( [&]( auto x ) { nodes.push_back( x ); } );
//   for (auto n : nodes) {
//     std::vector<uint64_t> fanin;
//     klut.foreach_fanin(n, [&](auto m) { fanin.push_back(m); });
//   fmt::print( "{} {}\n", n, fmt::join( fanin, " " ) );
//     for (auto m : fanin) {
//       fanout[m]++;
//     }
//   }

//  double splitter_cost = 2.0;
//  double buffer_cost = 2.0;

//   // computing aqfp cost
//   double extra_splitter_cost = 0.0;

//   for (auto n : nodes) {
//     if (klut.is_pi(n)) {
//       fmt::print("node {} is a pi\n", n);
//       continue;
//     }

//     if (klut.is_constant(n)) {
//       fmt::print("node {} is a constant\n", n);
//       continue;
//     }

//     if (fanout[n] > 1) {
//       extra_splitter_cost += num_splitters(fanout[n], 4u) * splitter_cost;
//     }

//     std::vector<uint64_t> fanin;
//     klut.foreach_fanin(n, [&](auto m) { fanin.push_back(m); });

//     if (fanin.size() == 1u) {
//       level[n] = level[fanin[0]] + 1;
//       continue;
//     }

//     std::vector<uint32_t> fanin_levels(4u, 0u);
//     for (auto i = 0u; i < fanin.size(); i++) {
//       fanin_levels[i] = level[fanin[i]] + (klut.is_constant(fanin[i]) ? 0u : num_splitter_levels(fanout[fanin[i]], 4u));
//     }

//     auto tt = kitty::extend_to(klut.node_function(n), 4u);
//     fmt::print("node {} to synthesize {:04x}\n", n, tt._bits[0]);
//     auto res = db.get_best_replacement(tt._bits[0], fanin_levels, buffer_cost);
//     fmt::print("can be synthesized at a cost of {}\n", std::get<2>(res));
//     fmt::print("fanin = {}\n", fmt::join(fanin, " "));
//     fmt::print("level configuration of inputs = {}\n", fmt::join(fanin_levels, " "));
//     fmt::print("level configuration of target = {}\n", fmt::join(std::get<1>(res).levels, " "));
//     fmt::print("best level = {}\n", std::get<3>(res));
//     total_cost += std::get<2>(res);
//     level[n] = std::get<3>(res);

//   }

//   total_cost += extra_splitter_cost;
//   fmt::print( "total cost = {}\n", total_cost );
// }

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;
  // mockturtle::aqfp_logical_network_db<mockturtle::aqfp_logical_network_t<int>, 4u> db( 2.0, 1u );
  //   std::ifstream dbchunk( "aqfp_db.txt" );
  //   assert( dbchunk.is_open() );
  //   db.load_db_from_file( dbchunk );
  // std::ofstream final_db( "aqfp_db_corrected.txt" );
  // db.save_db_to_file( final_db );
  // final_db.close();

  // test_partition_generation();
  // test_partition_extension();
  // test_sublist_generation();
  // test_layer_extension();
  // test_dag_simulation();
  // test_timing();
  // assemble_db();
  combine_db();
  // reverse_permutations();
  // test_aqfp_syn();

  return 0;
}
