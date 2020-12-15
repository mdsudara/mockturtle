#include <chrono>
#include <fmt/format.h>
#include <future>
#include <iostream>
#include <mutex>
#include <queue>
#include <sstream>
#include <thread>
#include <vector>

#include <mockturtle/algorithms/exact_syn/dag.hpp>
#include <mockturtle/algorithms/exact_syn/sat.hpp>

#include "./aqfp_cost_all_config.hpp"
#include "./aqfp_cost_func.hpp"
#include "./npn4.hpp"
#include "./simulate_dag.hpp"

template<typename T>
class task_queue
{
public:
  void push( T t )
  {
    m.lock();
    tasks.push( t );
    m.unlock();
  }

  T pop()
  {
    m.lock();
    auto res = tasks.front();
    tasks.pop();
    m.unlock();
    return res;
  }

  bool empty()
  {
    return tasks.empty();
  }

  size_t size()
  {
    return tasks.size();
  }

private:
  std::queue<T> tasks;
  std::mutex m;
};

using syn_task_t = std::tuple<mockturtle::aqfp_logical_network, size_t, uint64_t>;

task_queue<syn_task_t> syn_tasks;
std::atomic<size_t> num_syn_task_threads;

const size_t MAX_SYN_TASK_THREADS = 48;
const size_t MAX_SYN_TASK_QUEUE_SIZE = 100;

size_t num_in = 4u;

size_t to_synthesize = 0u;
std::atomic<size_t> num_synthesized = 0u;
std::vector<uint64_t> tt_to_synthesize = {};

std::vector<size_t> synthesized( ( 1u << 16 ), 0u );

/*
 * For each partial DAG, hold the list of functions that can be synthesized by that partial DAG.
 * For any DAG derived from a given partial DAG, we only need to check functions that could be 
 * synthesized with that partial DAG.
 */
std::vector<std::vector<uint64_t>> pdag_tts( 1000000u );

/**
 * Check if a dag can 'net' can synthesize given function 'tt' with 'num_in' inputs.
 */
void synthesize_async( const mockturtle::aqfp_logical_network& net, size_t num_in, uint64_t tt )
{
  assert( net.is_dag );
  mockturtle::sat_enc enc( net, num_in, tt, 0 );
  bool res = enc.solve();

  if ( res )
  {
    auto res = synthesized[tt]++;
    if ( !res )
    {
      num_synthesized++;
    }

    fmt::print( "{:04x} done with cost {:6.2f}\n{}\n", tt, net.computed_cost, enc.chain_as_string() );
  }

  num_syn_task_threads--;
}

/**
 * While there exists unprocessed synthesis tasks, process them in separate threads.
 */
void process_syn_tasks()
{
  fmt::print( "starting the process syn task thread\n" );

  while ( num_synthesized < to_synthesize )
  {
    if ( syn_tasks.empty() || num_syn_task_threads >= MAX_SYN_TASK_THREADS )
    {
      std::this_thread::sleep_for( std::chrono::milliseconds( 100 ) );
    }
    else
    {
      auto task = syn_tasks.pop();

      if ( synthesized[std::get<2>( task )] > 0 )
      {
        continue;
      }

      fmt::print( "synthesizing function {:04x} with dag {}\n", std::get<2>( task ), mockturtle::as_string( std::get<0>( task ) ) );

      num_syn_task_threads++;
      while ( true )
      {
        try
        {
          std::thread thrd( synthesize_async, std::get<0>( task ), std::get<1>( task ), std::get<2>( task ) );
          thrd.detach();
          break;
        }
        catch ( const std::exception& e )
        {
          fmt::print( "async error when synthesizing: {}\n", e.what() );
          std::this_thread::sleep_for( std::chrono::milliseconds( 10u ) );
        }
      }
    }
  }
}

void test_aqfp_costs()
{
  // size_t pdag_id = 0u;
  // bool is_partial_dag = 0u;
  // bool is_dag = 1u;
  // double computed_cost = -1.0;
  // size_t level = 0u;
  // std::vector<std::vector<node_t>> nodes;
  // std::vector<size_t> node_num_fanin;
  // std::map<size_t, size_t> num_gates_of_fanin;
  // std::vector<node_t> last_layer_leaves;
  // std::vector<node_t> other_leaves;
  // std::vector<node_t> input_slots;
  // size_t zero_input = 0;

  mockturtle::aqfp_logical_network net1 = {
      0u,
      false,
      true,
      -1.0,
      3u,
      { { 1, 4, 5 }, { 2, 4, 5 }, { 3, 5, 6 }, { 5, 7, 8 }, {}, {}, {}, {}, {} },
      { 3, 3, 3, 3, 0, 0, 0, 0, 0 },
      { { 3u, 4u } },
      {},
      {},
      { 4, 5, 6, 7, 8 },
      5u };

  fmt::print( "cost of net1 is {}\n", aqfp_cost( net1 ) );

  mockturtle::aqfp_logical_network net2 = {
      0u,
      false,
      true,
      -1.0,
      3u,
      { { 1, 2, 6 }, { 4, 5, 6 }, { 3, 6, 7 }, { 6, 8, 9 }, { 6, 7, 10 }, { 6, 8, 9 }, {}, {}, {}, {}, {} },
      { 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0 },
      { { 3u, 6u } },
      {},
      {},
      { 6, 7, 8, 9, 10 },
      0u };

  fmt::print( "cost of net2 is {}\n", aqfp_cost( net2 ) );
}

/**
 * Computes costs for DAGs in order. No synthesis is done.
 */
void test_dag_generation()
{
  std::vector<size_t> allowed_num_fanins = { 3u };                          // will use only fanin 3 gates
  std::map<size_t, size_t> max_gates_of_fanin = { { 3u, 7u }, { 5u, 0u } }; // allow at most 7 3-input gates and 0 5-input gates

  auto params = mockturtle::dag_params();

  params.max_gates = 7u;      // allow at most 7 gates in total
  params.max_num_fanout = 4u; // limit the maximum fanout of a gate
  params.max_width = 1000u;   // maximum number of gates at any level
  params.max_num_in = 5u;     // maximum number of input slots (need extra one for the constant)
  params.max_level = 4u;      // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  // auto gen = mockturtle::gen_dag( params, simple_cost );
  auto gen = mockturtle::gen_dag( params, aqfp_cost );

  auto count = 0u;
  while ( count++ < 1000000000u )
  {
    auto result_opt = gen.next_potentially_feasible( []( mockturtle::aqfp_logical_network result ) { (void) result; return true; } );

    if ( result_opt == std::nullopt )
    {
      break;
    }

    auto net = result_opt.value();

    fmt::print( "{} network {} has cost {}\n", count, mockturtle::as_string( net ), net.computed_cost );
  }
}

void save_all()
{
  std::vector<size_t> allowed_num_fanins = { 3u };              // will use only fanin 3 gates
  std::map<size_t, size_t> max_gates_of_fanin = { { 3u, 7u } }; // allow at most 7 3-input gates and 0 5-input gates

  auto params = mockturtle::dag_params();

  params.max_gates = 7u;         // allow at most 6 gates in total
  params.max_num_fanout = 1000u; // limit the maximum fanout of a gate
  params.max_width = 1000u;      // maximum number of gates at any level
  params.max_num_in = 5u;        // maximum number of inputs slots (need extra one for the constant)
  params.max_level = 7u;         // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  auto gen = mockturtle::gen_dag( params, simple_cost );

  auto t0 = std::chrono::high_resolution_clock::now();
  auto count = 0u;
  while ( count < 1000000000u )
  {
    auto result_opt = gen.next_potentially_feasible( []( auto& net ) {(void) net; return true; } );

    if ( result_opt == std::nullopt )
    {
      /* No more DAGs */
      break;
    }
    auto result = result_opt.value();

    fmt::print( "{}\n", result.get_encoding() );

    count++;
    if ( count % 1000 == 0 )
    {
      auto t1 = std::chrono::high_resolution_clock::now();
      auto d1 = std::chrono::duration_cast<std::chrono::microseconds>( t1 - t0 );
      std::cerr << fmt::format( "dags generated {:8d} time so far milliseconds {:8d}\n", count, d1.count() / 1000ul );
    }
  }
}

void synthesize_all()
{
  std::vector<size_t> allowed_num_fanins = { 3u };                          // will use only fanin 3 gates
  std::map<size_t, size_t> max_gates_of_fanin = { { 3u, 7u }, { 5u, 0u } }; // allow at most 7 3-input gates and 0 5-input gates

  auto params = mockturtle::dag_params();

  params.max_gates = 7u;         // allow at most 7 gates in total
  params.max_num_fanout = 1000u; // limit the maximum fanout of a gate
  params.max_width = 1000u;      // maximum number of gates at any level
  params.max_num_in = 5u;        // maximum number of inputs slots (need extra one for the constant)
  params.max_level = 4u;         // maximum number of gate levels in a DAG

  params.allowed_num_fanins = allowed_num_fanins;
  params.max_gates_of_fanin = max_gates_of_fanin;

  // auto gen = mockturtle::gen_dag( params, simple_cost );
  auto gen = mockturtle::gen_dag( params, aqfp_cost );

  auto count = 0u;
  while ( count < 1000000000u && num_synthesized < to_synthesize )
  {
    auto result_opt = gen.next_potentially_feasible(
        /* Return true if a partial DAG can synthesize any of the remaining functions. 
		   Otherwise, no need to consider the DAGs generated by it. */
        []( mockturtle::aqfp_logical_network result ) {
          fmt::print( "check feasibility of for {}\n", mockturtle::as_string( result ) );

          assert( result.is_partial_dag );
          for ( auto&& t : result.last_layer_leaves )
          {
            result.add_leaf_node( { t } );
          }

          for ( auto&& t : result.other_leaves )
          {
            result.add_leaf_node( { t } );
          }

          result.last_layer_leaves.clear();
          result.other_leaves.clear();

          std::vector<std::future<std::pair<bool, uint64_t>>> results;
          std::atomic<size_t> num_threads = 0u;

          for ( auto&& tt : tt_to_synthesize )
          {
            if ( synthesized[tt] > 0 )
            {
              continue;
            }

            while ( num_threads >= 48 )
            {
              std::this_thread::sleep_for( std::chrono::milliseconds( 1u ) );
            }

            num_threads++;

            while ( true )
            {
              try
              {
                results.push_back( std::async(
                    std::launch::async,
                    [&]( mockturtle::aqfp_logical_network& net, size_t num_in, uint64_t tt ) -> std::pair<bool, uint64_t> {
                      mockturtle::sat_enc enc( net, num_in, tt );
                      auto res = enc.solve();
                      if ( res )
                      {
                        // fmt::print( "pdag with id {} can synthesize {:04x}\n{}\n", net.pdag_id, tt, enc.chain_as_string() );
                      }
                      else
                      {
                        // fmt::print( "pdag cannot synthesize {:04x}\n", tt );
                      }
                      num_threads--;
                      return { res, tt };
                    },
                    std::ref( result ), num_in, tt ) );

                break;
              }
              catch ( const std::exception& e )
              {
                fmt::print( "async error when feasibility checking: {}\n", e.what() );
                std::this_thread::sleep_for( std::chrono::milliseconds( 10u ) );
              }
            }
          }

          bool ret = false;

          for ( auto&& t : results )
          {
            auto res = t.get();
            if ( res.first && synthesized[res.second] == 0u )
            {
              ret = true;
              pdag_tts[result.pdag_id].push_back( res.second );
            }
          }

          return ret;
        } );

    if ( result_opt == std::nullopt )
    {
      /* No more DAGs */
      break;
    }

    auto result = result_opt.value();

    if ( ( result.input_slots.size() > num_in + 1 ) || ( result.input_slots.size() == num_in + 1 && result.zero_input == 0u ) )
    {
      continue;
    }

    fmt::print( "{:8d} dag {} with cost {} potential funtions {}\n", ++count, mockturtle::as_string( result ), result.computed_cost, pdag_tts[result.pdag_id].size() );

    for ( auto&& tt : pdag_tts[result.pdag_id] )
    {
      while ( true )
      {
        if ( syn_tasks.size() < MAX_SYN_TASK_QUEUE_SIZE )
        {
          syn_tasks.push( { result, num_in, tt } );
          break;
        }
        else
        {
          std::this_thread::sleep_for( std::chrono::milliseconds( 100u ) );
        }
      }
    }

    /* Occasionally update the remaining functions to synthesize by removing the ones that have 
	   been synthesized so far. This will speed-up the feasibility checking for a partial DAG. */
    if ( count % 1000u == 0 )
    {
      std::vector<uint64_t> temp;
      for ( auto&& tt : tt_to_synthesize )
      {
        if ( synthesized[tt] == 0u )
        {
          temp.push_back( tt );
        }
      }
      tt_to_synthesize = temp;
      fmt::print( "updated tt-to-synthesize with new size {:3d}\n", tt_to_synthesize.size() );
    }

    if ( count % 100u == 0 )
    {
      fmt::print( "progress {:3d} out of {:3d}\n", num_synthesized, to_synthesize );
    }
  }
}

/**
 * Read an input stream until end and return the list of lines.
 */
std::vector<std::string> read_lines( std::istream& is )
{
  std::vector<std::string> result;
  std::string temp;
  while ( getline( is, temp ) )
  {
    if ( temp.length() > 0 )
    {
      result.push_back( temp );
    }
  }
  return result;
}

void process_all_dags()
{
  std::string temp;
  size_t count = 0u;
  while ( getline( std::cin, temp ) )
  {
    if ( temp.length() > 0 )
    {
      fmt::print( "processing dag {} [{}]\n", ++count, temp );
      mockturtle::aqfp_logical_network_t net;
      net.decode_string( temp );
      // fmt::print("{}\n", mockturtle::as_string(net));
      aqfp_cost_all( net );
    }
  }
}

int main( int argc, char** argv )
{
  (void)argc;
  (void)argv;

  // simulate_all_dags();
  // process_all_dags();
  save_all();

  return 0;

  for ( auto i = 0u; i < NUM_NPN_4; i++ )
  {
    if ( ( NPN_4[i] == 0x0000 ) || ( NPN_4[i] == 0x00ff ) )
    {
      continue;
    }
    tt_to_synthesize.push_back( NPN_4[i] );
  }

  to_synthesize = tt_to_synthesize.size();

  fmt::print( "number of functions to synthesize: {}\n", to_synthesize );

  /* Strat a new thread to process synthesis tasks. */
  std::thread syn_task_thread( process_syn_tasks );

  synthesize_all();

  /* Wait for syn_task_thread to finsih */
  syn_task_thread.join();

  while ( num_syn_task_threads > 0 )
  {
    fmt::print( "waiting for pending syn tasks to finish...\n" );
    std::this_thread::sleep_for( std::chrono::milliseconds( 100u ) );
  }

  return 0;
}
