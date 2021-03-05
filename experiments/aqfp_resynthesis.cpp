#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>
#include <lorina/lorina.hpp>

#include <mockturtle/io/blif_reader.hpp>
#include <mockturtle/io/verilog_reader.hpp>
#include <mockturtle/io/write_blif.hpp>

#include <mockturtle/networks/aqfp.hpp>
#include <mockturtle/networks/klut.hpp>
#include <mockturtle/networks/mig.hpp>

#include <mockturtle/views/depth_view.hpp>
#include <mockturtle/views/fanout_limit_view.hpp>
#include <mockturtle/views/fanout_view.hpp>

#include <mockturtle/algorithms/aqfp_resynthesis.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_db_builder.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_fanout_resyn.hpp>
#include <mockturtle/algorithms/aqfp_resynthesis/aqfp_node_resyn.hpp>

#include <mockturtle/properties/aqfpcost.hpp>

#include "../experiments/experiments.hpp"

template<typename Result>
bool has_better_cost( Result& current, Result& previous )
{
  if ( current.first < previous.first )
    return true;

  if ( current.first > previous.first )
    return false;

  return current.second < previous.second;
}

template<typename Result>
bool has_better_level( Result& current, Result& previous )
{
  if ( current.second < previous.second )
    return true;

  if ( current.second > previous.second )
    return false;

  return current.first < previous.first;
}

std::vector<std::string> mcnc = {
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
mockturtle::klut_network lut_map( Ntk const& ntk, uint32_t k = 4 )
{
  std::string tempfile1 = "/tmp/temp1.blif";
  std::string tempfile2 = "/tmp/temp2.blif";

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

void experiment_aqfp_exact_syn( std::unordered_map<uint32_t, double>& gate_costs, std::unordered_map<uint32_t, double>& splitters, mockturtle::aqfp_db<>& db, const std::vector<std::string>& benchmarks )
{
  mockturtle::aqfp_network_cost cost_fn( gate_costs, splitters );

  mockturtle::aqfp_node_resyn node_resyn_cst( db, { splitters, mockturtle::aqfp_node_resyn_strategy::cost_based } );
  mockturtle::aqfp_node_resyn node_resyn_lvl( db, { splitters, mockturtle::aqfp_node_resyn_strategy::level_based } );
  mockturtle::aqfp_fanout_resyn fanout_resyn( 4u );

  experiments::experiment<std::string, double, uint32_t, double, uint32_t, double, uint32_t, double, uint32_t> exp( "aqfp_resynthesis", "benchmark", "#JJ (C01)", "LVL (C01)", "#JJ (C10)", "LVL (C10)", "#JJ (L01)", "LVL (L01)", "#JJ (L10)", "LVL (L10)" );

  for ( auto b : benchmarks )
  {
    fmt::print( "Processing benchmark {}...\n", b );

    std::string benchmark = fmt::format( "./benchmarks/{}", b );
    mockturtle::mig_network mig;

    lorina::read_verilog( benchmark, mockturtle::verilog_reader( mig ) );
    fmt::print( "\tpi: {:4d} po: {:4d} size: {:6d}\n", mig.num_pis(), mig.num_pos(), mig.num_gates() );

    mockturtle::klut_network klut_orig = lut_map( mig, 4 );

    /* 1. Apply cost-based AQFP resynthesis once */
    mockturtle::aqfp_network opt_aqfp;
    auto res = mockturtle::aqfp_resynthesis( opt_aqfp, klut_orig, node_resyn_cst, fanout_resyn );
    std::pair<double, uint32_t> res_orig_cst = { cost_fn( opt_aqfp, res.node_level, res.po_level ), res.po_level };

    /* 2. Repeatedly apply cost-based AQFP resynthesis */
    auto res_opt_cst = res_orig_cst;
    for ( auto i = 2u; i <= 10u; i++ )
    {
      auto klut_opt = lut_map( opt_aqfp, 4 );

      opt_aqfp = mockturtle::aqfp_network();
      res = mockturtle::aqfp_resynthesis( opt_aqfp, klut_opt, node_resyn_cst, fanout_resyn );
      std::pair<double, uint32_t> res_temp = { cost_fn( opt_aqfp, res.node_level, res.po_level ), res.po_level };

      if ( has_better_cost( res_temp, res_opt_cst ) )
      {
        res_opt_cst = res_temp;
      }
    }

    assert( abc_cec_aqfp( opt_aqfp, benchmark ) );

    /* 3. Apply level-based AQFP resynthesis once */
    opt_aqfp = mockturtle::aqfp_network();
    res = mockturtle::aqfp_resynthesis( opt_aqfp, klut_orig, node_resyn_lvl, fanout_resyn );
    std::pair<double, uint32_t> res_orig_lvl = { cost_fn( opt_aqfp, res.node_level, res.po_level ), res.po_level };

    /* 4. Repeatedly apply level-based AQFP resynthesis */
    auto res_opt_lvl = res_orig_lvl;
    for ( auto i = 2u; i <= 10u; i++ )
    {
      auto klut_opt = lut_map( opt_aqfp, 4 );

      opt_aqfp = mockturtle::aqfp_network();
      res = mockturtle::aqfp_resynthesis( opt_aqfp, klut_opt, node_resyn_lvl, fanout_resyn );
      std::pair<double, uint32_t> res_temp = { cost_fn( opt_aqfp, res.node_level, res.po_level ), res.po_level };

      if ( has_better_level( res_temp, res_opt_lvl ) )
      {
        res_opt_lvl = res_temp;
      }
    }

    assert( abc_cec_aqfp( opt_aqfp, benchmark ) );

    exp( b, res_orig_cst.first, res_orig_cst.second, res_opt_cst.first, res_opt_cst.second, res_orig_lvl.first, res_orig_lvl.second, res_opt_lvl.first, res_opt_lvl.second );
  }

  exp.save();
  exp.table();
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  std::unordered_map<uint32_t, double> gate_costs = { { 3u, 6.0 }, { 5u, 10.0 } };
  std::unordered_map<uint32_t, double> splitters = { { 1u, 2.0 }, { 4u, 2.0 } };

  mockturtle::aqfp_db_builder<> builder( gate_costs, splitters );

  std::ifstream db_file( "aqfp_db.txt" );
  assert( db_file.is_open() );
  builder.load_db_from_file( db_file );
  db_file.close();

  // mockturtle::aqfp_db<> db( gate_costs, splitters );
  auto db = builder.build();

  experiment_aqfp_exact_syn( gate_costs, splitters, db, mcnc );

  return 0;
}
