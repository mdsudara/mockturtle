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

#include <mockturtle/algorithms/exact_syn/aqfp_syn.hpp>
#include <mockturtle/algorithms/exact_syn/dag.hpp>

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

bool has_better_cost( mockturtle::aqfp_syn_result& current, mockturtle::aqfp_syn_result& previous )
{
  if ( current.total_cost < previous.total_cost )
    return true;

  if ( current.total_cost > previous.total_cost )
    return false;

  return current.max_level < previous.max_level;
}

bool has_better_level( mockturtle::aqfp_syn_result& current, mockturtle::aqfp_syn_result& previous )
{
  if ( current.max_level < previous.max_level )
    return true;

  if ( current.max_level > previous.max_level )
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

    "chkn.v",
    "count.v",
    "dist.v",
    "in5.v",
    "in6.v",
    "k2.v",
    "m3.v",
    "max512.v",
    "misex3.v",
    "mlp4.v",
    "prom2.v",
    "sqr6.v",
    "x1dn.v",
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
void experiment_aqfp_exact_syn_only( AqfpSynT& aqfp_syn, const std::vector<std::string>& benchmarks, const std::string& experiment_name, const std::string& strategy )
{
  experiments::experiment<
      std::string,
      uint32_t, double,
      uint32_t, double,
      uint32_t, double,
      uint32_t, double>
      exp(
          experiment_name,
          "benchmark",
          "    JJ Depth", "        #JJ",
          " FL JJ Depth", "     FL #JJ",
          "OPT JJ Depth", "    OPT #JJ",
          "OPT FL Depth", "OPTL FL #JJ" );

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
    auto res_orig = aqfp_syn.synthesize( klut_orig, aqfp_net_1 );

    /* 2. Limit fan-outs and then map to AQFP */
    mockturtle::mig_network fo_lim_mig;
    mockturtle::fanout_limit_view fo_lim_v( fo_lim_mig, mockturtle::fanout_limit_view_params{ 16u } );
    cleanup_dangling( mig, fo_lim_v );

    mockturtle::klut_network klut_folim = lut_map_old( fo_lim_v, benchmark, "temp", 4 );
    mockturtle::aqfp_network aqfp_net_2;
    auto res_folim = aqfp_syn.synthesize( klut_folim, aqfp_net_2 );

    /* 3. Repeatedly apply exact AQFP synthesis */
    mockturtle::aqfp_network opt_aqfp;
    auto res_opt = aqfp_syn.synthesize( klut_orig, opt_aqfp );

    mockturtle::aqfp_network opt_folim_aqfp;
    auto res_folim_opt = aqfp_syn.synthesize( klut_folim, opt_folim_aqfp );

    for ( auto i = 2u; i <= 10u; i++ )
    {
      auto klut_opt = lut_map_old( opt_aqfp, benchmark, "temp", 4 );

      mockturtle::mig_network temp_mig;
      mockturtle::fanout_limit_view temp_folim_mig( temp_mig, mockturtle::fanout_limit_view_params{ 16u } );
      cleanup_dangling( opt_folim_aqfp, temp_folim_mig );

      auto klut_folim_opt = lut_map_old( temp_folim_mig, benchmark, "temp", 4 );

      opt_aqfp = mockturtle::aqfp_network();
      auto res_opt1 = aqfp_syn.synthesize( klut_opt, opt_aqfp );

      opt_folim_aqfp = mockturtle::aqfp_network();
      auto res_folim_opt1 = aqfp_syn.synthesize( klut_folim_opt, opt_folim_aqfp );

      if ( strategy == "cost" && has_better_cost( res_opt1, res_opt ) )
      {
        res_opt = res_opt1;
      }
      if ( strategy == "level" && has_better_level( res_opt1, res_opt ) )
      {
        res_opt = res_opt1;
      }

      if ( strategy == "cost" && has_better_cost( res_folim_opt1, res_folim_opt ) )
      {
        res_folim_opt = res_folim_opt1;
      }
      if ( strategy == "level" && has_better_cost( res_folim_opt1, res_folim_opt ) )
      {
        res_folim_opt = res_folim_opt1;
      }
    }

    assert( abc_cec_aqfp( opt_aqfp, benchmark ) );
    assert( abc_cec_aqfp( opt_folim_aqfp, benchmark ) );

    exp( b,
         res_orig.max_level, res_orig.total_cost,
         res_folim.max_level, res_folim.total_cost,
         res_opt.max_level, res_opt.total_cost,
         res_folim_opt.max_level, res_folim_opt.total_cost );
  }

  exp.save( experiment_name );

  fmt::print( "Experiment details\n"
              "\tAQFP exact-syn on original network\n"
              "\tAQFP exact-syn on fanout-limited original network\n"
              "\tAQFP exact-syn best out of 10 iterations\n"
              "\tAQFP exact-syn best out of 10 iterations with fanout limitting\n" );
  exp.table();
}

template<typename SrcNtk>
mockturtle::mig_network selective_depth_rewrite( SrcNtk& ntk )
{
  mockturtle::depth_view_params dv_params;
  mockturtle::depth_view dv{ ntk, fanout_cost_depth_local<mockturtle::mig_network>(), dv_params };

  mockturtle::mig_algebraic_depth_rewriting_params ps_a;
  ps_a.strategy = mockturtle::mig_algebraic_depth_rewriting_params::selective;
  ps_a.allow_area_increase = true;

  mig_algebraic_depth_rewriting_splitters( dv, ps_a );

  mockturtle::mig_network mig;
  cleanup_dangling( ntk, mig );

  return mig;
}

template<typename SrcNtk>
mockturtle::mig_network splitter_aware_resubstitute( SrcNtk& ntk )
{
  mockturtle::resubstitution_params ps_r;
  ps_r.max_divisors = 250;
  ps_r.max_inserts = 1;

  mockturtle::fanout_view fv( ntk );

  mockturtle::depth_view_params dv_params;
  mockturtle::depth_view dv( fv, fanout_cost_depth_local<mockturtle::mig_network>(), dv_params );

  mig_resubstitution_splitters( dv, ps_r );

  mockturtle::mig_network mig;
  cleanup_dangling( dv, mig );

  return mig;
}

template<typename SrcNtk>
mockturtle::mig_network splitter_aware_ackers_resyn( SrcNtk& ntk )
{
  //  mockturtle::fanout_view fv( ntk );

  mockturtle::akers_resynthesis<mockturtle::mig_network> resyn;
  refactoring( ntk, resyn, {}, nullptr, jj_cost<mockturtle::mig_network>() );

  mockturtle::mig_network mig;
  cleanup_dangling( ntk, mig );

  return mig;
}

template<typename SrcNtk>
mockturtle::mig_network splitter_aware_npn_resyn( SrcNtk& ntk )
{
  //  mockturtle::fanout_view fv( ntk );

  mockturtle::mig_npn_resynthesis resyn;

  mockturtle::refactoring_params ps_r;
  ps_r.max_pis = 4u;

  refactoring( ntk, resyn, ps_r, nullptr, jj_cost<mockturtle::mig_network>() );

  mockturtle::mig_network mig;
  cleanup_dangling( ntk, mig );

  return mig;
} 

template<typename AqfpSynT>
void experiment_aqfp_exact_syn_with_asp_dac_flow( AqfpSynT& aqfp_syn, const std::vector<std::string>& benchmarks, const std::string& experiment_name, const std::string& strategy )
{
  experiments::experiment<std::string, uint32_t, double, uint32_t, double> exp(
      experiment_name,
      "benchmark",
      "    JJ Depth", "        #JJ",
      "OPT JJ Depth", "    OPT #JJ" );

  for ( auto b : benchmarks )
  {
    fmt::print( "Processing benchmark {}...\n", b );

    std::string benchmark = fmt::format( "./benchmarks/{}", b );
    mockturtle::mig_network mig;

    lorina::read_verilog( benchmark, mockturtle::verilog_reader( mig ) );
    fmt::print( "\tpi: {:4d} po: {:4d} size: {:6d}\n", mig.num_pis(), mig.num_pos(), mig.num_gates() );

    /* 1. Directly map to AQFP */
    mockturtle::klut_network klut_orig = lut_map_old( mig, benchmark, "temp", 4 );
    mockturtle::aqfp_network opt_aqfp;
    auto res_orig = aqfp_syn.synthesize( klut_orig, opt_aqfp );

    mockturtle::aqfp_syn_result res_opt{
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<uint32_t>::max() };

    for ( auto i = 0u; i < 5u; i++ )
    {
      // convert aqfp network to an mig

      mockturtle::mig_network temp_mig;
       fmt::print("cleanup start\n");
      cleanup_dangling( opt_aqfp, temp_mig );

       fmt::print("rewrite start\n");
      temp_mig = selective_depth_rewrite( temp_mig );

       fmt::print("resub start\n");
      temp_mig = splitter_aware_resubstitute( temp_mig );
       fmt::print("resyn start\n");
      temp_mig = splitter_aware_npn_resyn( temp_mig );
//       temp_mig = splitter_aware_ackers_resyn( temp_mig );
       fmt::print("opt end\n");

      mig = temp_mig;

      mockturtle::klut_network klut_temp = lut_map_old( temp_mig, benchmark, "temp", 4 );
      opt_aqfp = mockturtle::aqfp_network();
      auto res_temp = aqfp_syn.synthesize( klut_temp, opt_aqfp );

      assert( abc_cec_aqfp( opt_aqfp, benchmark ) );

      fmt::print( "depth = {} cost = {}\n", res_temp.max_level, res_temp.total_cost );

      if ( has_better_cost(res_temp, res_opt) )
      {
        res_opt = res_temp;
      }
    }

    // klut_orig = lut_map_old( mig, benchmark, "temp", 4 );

    // opt_aqfp = mockturtle::aqfp_network();
    // res_opt = aqfp_syn.synthesize( klut_orig, opt_aqfp );

    // for ( auto i = 2u; i <= 10u; i++ )
    // {
    //   auto klut_opt = lut_map_old( opt_aqfp, benchmark, "temp", 4 );

    //   opt_aqfp = mockturtle::aqfp_network();
    //   auto res_opt1 = aqfp_syn.synthesize( klut_opt, opt_aqfp );

    //   if ( strategy == "cost" && has_better_cost( res_opt1, res_opt ) )
    //   {
    //     res_opt = res_opt1;
    //   }
    //   if ( strategy == "level" && has_better_level( res_opt1, res_opt ) )
    //   {
    //     res_opt = res_opt1;
    //   }
    // }

    // /* 3. Repeatedly apply exact AQFP synthesis */
    // mockturtle::aqfp_network opt_aqfp;
    // auto res_opt = aqfp_syn.synthesize( klut_orig, opt_aqfp );

    // mockturtle::aqfp_network opt_folim_aqfp;
    // auto res_folim_opt = aqfp_syn.synthesize( klut_folim, opt_folim_aqfp );

    // for ( auto i = 2u; i <= 10u; i++ )
    // {
    //   auto klut_opt = lut_map_old( opt_aqfp, benchmark, "temp", 4 );

    //   mockturtle::mig_network temp_mig;
    //   mockturtle::fanout_limit_view temp_folim_mig( temp_mig, mockturtle::fanout_limit_view_params{ 16u } );
    //   cleanup_dangling( opt_folim_aqfp, temp_folim_mig);

    //   auto klut_folim_opt = lut_map_old( temp_folim_mig, benchmark, "temp", 4 );

    //   opt_aqfp = mockturtle::aqfp_network();
    //   auto res_opt1 = aqfp_syn.synthesize( klut_opt, opt_aqfp );

    //   opt_folim_aqfp = mockturtle::aqfp_network();
    //   auto res_folim_opt1 = aqfp_syn.synthesize( klut_folim_opt, opt_folim_aqfp);

    //   if ( strategy == "cost" && ( res_opt1.total_cost < res_opt.total_cost || ( res_opt1.total_cost == res_opt.total_cost && res_opt1.max_level < res_opt.max_level ) ) )
    //   {
    //     res_opt = res_opt1;
    //   }
    //   if ( strategy == "level" && ( res_opt1.max_level < res_opt.max_level || ( res_opt1.max_level == res_opt.max_level && res_opt1.total_cost < res_opt.total_cost ) ) )
    //   {
    //     res_opt = res_opt1;
    //   }

    //   if ( strategy == "cost" && ( res_folim_opt1.total_cost < res_folim_opt.total_cost || ( res_folim_opt1.total_cost == res_folim_opt.total_cost && res_folim_opt1.max_level < res_folim_opt.max_level ) ) )
    //   {
    //     res_folim_opt = res_folim_opt1;
    //   }
    //   if ( strategy == "level" && ( res_folim_opt1.max_level < res_folim_opt.max_level || ( res_folim_opt1.max_level == res_folim_opt.max_level && res_folim_opt1.total_cost < res_folim_opt.total_cost ) ) )
    //   {
    //     res_folim_opt = res_folim_opt1;
    //   }
    // }

    // assert( abc_cec_aqfp( opt_aqfp, benchmark ) );
    // assert( abc_cec_aqfp( opt_folim_aqfp, benchmark ) );

    exp( b,
         res_orig.max_level, res_orig.total_cost,
         res_opt.max_level, res_opt.total_cost );
  }

  exp.save( experiment_name );

  fmt::print( "Experiment details\n"
              "\tAQFP exact-syn on original network\n"
              "\tAQFP exact-syn with algebraic rewriting - best of 10 iterations\n" );
  exp.table();
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  std::string method( argv[1] );

  mockturtle::aqfp_node_resyn_param ps_cost{ "aqfp_db_2020_feb_05.txt", { { 4u, 2.0 } }, 2.0, mockturtle::aqfp_node_resyn_strategy::best_cost, 0u };
  mockturtle::aqfp_node_resyn<mockturtle::aqfp_logical_network_t<int>> aqfp_syn_cost( ps_cost );

  if ( method == "cost" )
  {
    experiment_aqfp_exact_syn_only( aqfp_syn_cost, mcnc, "cost-based-aqfp-exact-syn-only-06-02-2021", method );

   // experiment_aqfp_exact_syn_with_asp_dac_flow( aqfp_syn_cost, mcnc, "cost-based-aqfp-exact-syn-with-asp-dac-flow-06-02-2021", method );
  }

   mockturtle::aqfp_node_resyn_param ps_level{ "aqfp_db_2020_feb_05.txt", { { 4u, 2.0 } }, 2.0, mockturtle::aqfp_node_resyn_strategy::best_level, 0u };
   mockturtle::aqfp_node_resyn<mockturtle::aqfp_logical_network_t<int>> aqfp_syn_level( ps_level );

   if ( method == "level" )
   {
      experiment_aqfp_exact_syn_only( aqfp_syn_level, mcnc, "level-based-aqfp-opt-06-02-2021", method );
   // experiment_aqfp_exact_syn_with_asp_dac_flow( aqfp_syn_level, mcnc, "level-based-aqfp-opt-06-02-2021", method );
   }

  return 0;
}
