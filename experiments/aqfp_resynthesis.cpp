#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>

#include <lorina/lorina.hpp>
#include <mockturtle/generators/majority.hpp>
#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>
#include <mockturtle/networks/aqfp.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_limit_view.hpp>
#include <mockturtle/views/fanout_view.hpp>

#include <mockturtle/algorithms/cleanup.hpp>

#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_node_resyn.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_fanout_resyn.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/dag.hpp>

#include <mockturtle/algorithms/aqfp_resynthesis.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting.hpp>
#include <mockturtle/algorithms/mig_algebraic_rewriting_splitters.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/mig_resub_splitters.hpp>
#include <mockturtle/algorithms/node_resynthesis.hpp>
#include <mockturtle/algorithms/node_resynthesis/akers.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/mig_npn.hpp>
#include <mockturtle/algorithms/refactoring.hpp>
#include <mockturtle/views/fanout_limit_view.hpp>

#include "../experiments/experiments.hpp"

template<typename Result>
bool has_better_cost( Result& current, Result& previous )
{
  if ( current.total_cost < previous.total_cost )
    return true;

  if ( current.total_cost > previous.total_cost )
    return false;

  return current.critical_po_level < previous.critical_po_level;
}

template<typename Result>
bool has_better_level( Result& current, Result& previous )
{
  if ( current.critical_po_level < previous.critical_po_level )
    return true;

  if ( current.critical_po_level > previous.critical_po_level )
    return false;

  return current.total_cost < previous.total_cost;
}

template<class Ntk>
struct jj_cost
{
  uint32_t operator()( Ntk const& ntk, mockturtle::node<Ntk> const& n ) const
  {
    if ( ntk.is_pi( n ) || ntk.is_constant( n ) )
      return 0;

    return 6 + 2 * ( ( ntk.fanout_size( n ) + 1 ) / 3 );
  }
};

template<class Ntk>
struct fanout_cost_depth_local
{
  uint32_t operator()( Ntk const& ntk, mockturtle::node<Ntk> const& n ) const
  {
    uint32_t cost = 0;
    if ( ntk.is_pi( n ) )
      cost = 0;
    else if ( ntk.fanout_size( n ) == 1 )
      cost = 1;
    else if ( ( ntk.fanout_size( n ) > 1 ) && ( ( ntk.fanout_size( n ) < 5 ) ) )
      cost = 2;
    else if ( ( ntk.fanout_size( n ) > 4 ) && ( ntk.fanout_size( n ) < 17 ) )
      cost = 3;
    else
      cost = 4;
    return cost;
  }
};

std::vector<std::string> mcnc = {
    // "9symml.v",
    // "Z9sym.v",
    // "alu4.v",
    // "apex6.v",
    // "c17.v",
    // "c2670.v",
    // "c3540.v",
    // "c7552.v",
    // "ex7.v",
    // "in7.v",
    // "x4.v",
    // "x9dn.v",
    // "xparc.v",

    "5xp1.v",
    "c1908.v",
    "c432.v",
    "c5315.v",
    "c880.v",

    // "chkn.v",
    // "count.v",
    // "dist.v",
    // "in5.v",
    // "in6.v",
    // "k2.v",
    // "m3.v",
    // "max512.v",
    // "misex3.v",
    // "mlp4.v",
    // "prom2.v",
    // "sqr6.v",
    // "x1dn.v",
};

template<class Ntk>
bool abc_cec_aqfp( Ntk const& ntk, std::string const& benchmark )
{
  mockturtle::write_bench( ntk, "/tmp/test.bench" );
  std::string command = fmt::format( "abc -q \"cec -n {} /tmp/test.bench\"", benchmark );

  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, decltype( &pclose )> pipe( popen( command.c_str(), "r" ), pclose );
  if ( !pipe )
  {
    throw std::runtime_error( "popen() failed" );
  }
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    result += buffer.data();
  }

  if ( result.size() >= 23 && result.substr( 0u, 23u ) == "Networks are equivalent" )
  {
    fmt::print( "SUCCESS\n" );
    return true;
  }
  return false;
}

template<typename Ntk>
mockturtle::klut_network lut_map_old( Ntk const& ntk, std::string design, std::string command, uint32_t k = 4 )
{
  std::string tempfile1 = design + ".temp1old." + command + ".blif";
  std::string tempfile2 = design + ".temp2old." + command + ".blif";

  mockturtle::write_blif( ntk, tempfile1 );
  system( fmt::format( "abc -q \"{}; if -K {}; write_blif {}\" >> /dev/null 2>&1", tempfile1, k, tempfile2 ).c_str() );

  mockturtle::klut_network klut;
  if ( lorina::read_blif( tempfile2, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "FATAL OLD LUT MAP - Reading mapped network failed!" << std::endl;
    std::abort();
    return klut;
  }

  system( fmt::format( "rm {}", tempfile1 ).c_str() );
  system( fmt::format( "rm {}", tempfile2 ).c_str() );
  return klut;
}

template<typename Ntk>
mockturtle::klut_network lut_map_new( Ntk const& ntk, std::string design, std::string command, uint32_t k = 4 )
{
  std::string tempfile1 = design + ".temp1new." + command + ".blif";
  std::string tempfile2 = design + ".temp2new." + command + ".blif";

  mockturtle::write_blif( ntk, tempfile1 );
  system( fmt::format( "abc -q \"{}; &get; &if -K {}; &put; write_blif {}\" >> /dev/null 2>&1 ", tempfile1, k, tempfile2 ).c_str() );

  mockturtle::klut_network klut;
  if ( lorina::read_blif( tempfile2, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
  {
    std::cout << "FATAL NEW LUT MAP - Reading mapped network failed!" << std::endl;
    std::abort();
    return klut;
  }

  system( fmt::format( "rm {}", tempfile1 ).c_str() );
  system( fmt::format( "rm {}", tempfile2 ).c_str() );
  return klut;
}

template<typename AqfpSynT>
void experiment_aqfp_exact_syn_only( AqfpSynT& aqfp_resyn, const std::vector<std::string>& benchmarks, const std::string& experiment_name, const std::string& strategy )
{
  experiments::experiment<
      std::string,
      uint32_t, double,
      uint32_t, double>
      exp(
          experiment_name,
          "benchmark",
          "    JJ Depth", "        #JJ",
          "OPT JJ Depth", "    OPT #JJ" );

   mockturtle::aqfp_fanout_resyn fanout_resyn(4u);

  for ( auto b : benchmarks )
  {
    fmt::print( "Processing benchmark {}...\n", b );

    std::string benchmark = fmt::format( "./benchmarks/{}", b );
    mockturtle::mig_network mig;

    lorina::read_verilog( benchmark, mockturtle::verilog_reader( mig ) );
    fmt::print( "\tpi: {:4d} po: {:4d} size: {:6d}\n", mig.num_pis(), mig.num_pos(), mig.num_gates() );

    /* 1. Directly map to AQFP */
    mockturtle::klut_network klut_orig = lut_map_old( mig, benchmark, "temp", 4 );
    mockturtle::aqfp_network aqfp_net_1;
    auto res_orig = mockturtle::aqfp_resynthesis( aqfp_net_1, klut_orig, aqfp_resyn, fanout_resyn);

    /* 2. Repeatedly apply exact AQFP synthesis */
    mockturtle::aqfp_network opt_aqfp;
    auto res_opt = mockturtle::aqfp_resynthesis( opt_aqfp, klut_orig, aqfp_resyn, fanout_resyn );

    for ( auto i = 2u; i <= 10u; i++ )
    {
      auto klut_opt = lut_map_old( opt_aqfp, benchmark, "temp", 4 );

      opt_aqfp = mockturtle::aqfp_network();
      auto res_opt1 = mockturtle::aqfp_resynthesis( opt_aqfp, klut_opt, aqfp_resyn, fanout_resyn );

      if ( strategy == "cost" && has_better_cost( res_opt1, res_opt ) )
      {
        res_opt = res_opt1;
      }
      if ( strategy == "level" && has_better_level( res_opt1, res_opt ) )
      {
        res_opt = res_opt1;
      }
    }


    assert( abc_cec_aqfp( opt_aqfp, benchmark ) );

    exp( b,
         res_orig.critical_po_level, res_orig.total_cost,
         res_opt.critical_po_level, res_opt.total_cost );
  }

  exp.save( experiment_name );

  fmt::print( "Experiment details\n"
              "\tAQFP exact-syn on original network\n"
              "\tAQFP exact-syn best out of 10 iterations\n" );
  exp.table();
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  std::string method( argv[1] );

  mockturtle::aqfp_node_resyn_param ps_cost{ "aqfp_db_2020_feb_05.txt", { { 4u, 2.0 } }, 2.0, mockturtle::aqfp_node_resyn_strategy::best_cost, 0u };
  mockturtle::aqfp_node_resyn<mockturtle::node<mockturtle::aqfp_network>, mockturtle::node<mockturtle::klut_network>> aqfp_syn_cost( ps_cost );

  if ( method == "cost" )
  {
    experiment_aqfp_exact_syn_only( aqfp_syn_cost, mcnc, "cost-based-aqfp", method );
  }

  mockturtle::aqfp_node_resyn_param ps_level{ "aqfp_db_2020_feb_05.txt", { { 4u, 2.0 } }, 2.0, mockturtle::aqfp_node_resyn_strategy::best_level, 0u };
  mockturtle::aqfp_node_resyn<mockturtle::node<mockturtle::aqfp_network>, mockturtle::node<mockturtle::klut_network>> aqfp_syn_level( ps_level );

  if ( method == "level" )
  {
    experiment_aqfp_exact_syn_only( aqfp_syn_level, mcnc, "level-based-aqfp", method );
  }

  return 0;
}
