#include <string>
#include <vector>

#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/algorithms/mig_resub.hpp>
#include <mockturtle/algorithms/cleanup.hpp>
#include <mockturtle/algorithms/node_resynthesis/exact.hpp>
#include <mockturtle/algorithms/node_resynthesis/tig_exact.hpp>
#include <mockturtle/io/aiger_reader.hpp>
#include <mockturtle/networks/mig.hpp>
#include <mockturtle/networks/tig_dot.hpp>
#include <mockturtle/networks/tig_maj.hpp>
#include <mockturtle/networks/tig_onehot.hpp>
#include <mockturtle/networks/tig_gamble.hpp>
#include <mockturtle/networks/tig_and3.hpp>
#include <mockturtle/networks/tig_xorand.hpp>
#include <mockturtle/networks/tig_andxor.hpp>
#include <mockturtle/networks/tig_orand.hpp>
#include <mockturtle/networks/tig_mux.hpp>




#include <mockturtle/mockturtle.hpp>
#include "utils.hpp"

using namespace mockturtle;
using namespace std;

template <typename Ntk>
void print_before(Ntk&& ntk) {
	fmt::print("before        : pis = {:6d}  pos = {:6d}  nodes = {:7d}  levels = {:6d}\n", ntk.num_pis(), ntk.num_pos(), ntk.num_gates(), depth(ntk));
}

template <typename Ntk>
void print_lutmap(Ntk&& ntk, int iter) {
	fmt::print("lutmap iter {:2d}: pis = {:6d}  pos = {:6d}  nodes = {:7d}  levels = {:6d}\n", iter, ntk.num_pis(), ntk.num_pos(), ntk.num_gates(), depth(ntk));
}

template <typename Ntk>
void print_after(Ntk&& ntk, int iter) {
	fmt::print("after  iter {:2d}: pis = {:6d}  pos = {:6d}  nodes = {:7d}  levels = {:6d}\n", iter, ntk.num_pis(), ntk.num_pos(), ntk.num_gates(), depth(ntk));
}

/* runs a cut-rewrite iteration on ntk, stops if size does not decrease (for verilog)*/
template <typename Ntk, typename Resyn>
void cut_rewrite_iter_verilog(Ntk ntk, std::string path, std::string outpath, Resyn&& resyn, uint32_t maxiter = 1000) {
	uint64_t prev_size = 1000000000u; 
	uint64_t prev_levels = 1000000000u;

	cut_rewriting_params crp;
	crp.cut_enumeration_ps.cut_size = 4;	

	mockturtle::write_blif(ntk, outpath);
	for (auto i = 1u; i <= maxiter; i++) {
		if (i == 1u) {
			lorina::read_verilog(path, mockturtle::verilog_reader(ntk));
			print_before(ntk);

			prev_size = ntk.num_gates();
			prev_levels = depth(ntk);
			mockturtle::write_blif(ntk, outpath);
		}

		ntk = cut_rewriting(ntk, resyn, crp);
		ntk = cleanup_dangling(ntk);
		print_after(ntk, i + 1);

		if (ntk.num_gates() >= prev_size) {
			fmt::print("optimized network stats gates {:7d} levels {:7d}\n", prev_size, prev_levels);
			klut_network mapped;
		    mapped= save_mapped_old(outpath, outpath.substr(0, outpath.length() - 5) + "_lut6_old.blif");
			fmt::print("mappedold network stats gates {:7d} levels {:7d}\n", mapped.size(), depth(mapped));
		    mapped= save_mapped_new(outpath, outpath.substr(0, outpath.length() - 5) + "_lut6_new.blif");
			fmt::print("mappednew network stats gates {:7d} levels {:7d}\n", mapped.size(), depth(mapped));
			return;
		}

		prev_size = ntk.num_gates();
		prev_levels = depth(ntk);
		mockturtle::write_blif(ntk, outpath);
	}
}

/* runs a node resynthesis iteration on ntk, stops if size does not decrease (for verilog)*/
template <typename Ntk, typename Resyn, typename MapFn, typename SaveFn>
void resyn_iter_verilog(Ntk ntk, std::string path, std::string outpath, std::string name, Resyn&& resyn, MapFn&& mapfn, SaveFn&& savefn, uint32_t maxiter = 1000) {
	uint64_t prev_size = 1000000000u; 
	uint64_t prev_levels = 1000000000u;

	for (auto i = 1u; i <= maxiter; i++) {
		if (i == 1u) {
			lorina::read_verilog(path, mockturtle::verilog_reader(ntk));
			print_before(ntk);

			prev_size = ntk.num_gates();
			prev_levels = depth(ntk);
			mockturtle::write_blif(ntk, outpath);
		}

		klut_network klut = mapfn(ntk, path, name, 4); 	
		print_lutmap(klut, i);

		ntk = node_resynthesis<Ntk>(klut, resyn);
		ntk = cleanup_dangling(ntk);
		print_after(ntk, i);
		
		if (ntk.num_gates() >= prev_size) {
			fmt::print("optimized network stats gates {:7d} levels {:7d}\n", prev_size, prev_levels);
			klut_network mapped = savefn(outpath, outpath.substr(0, outpath.length() - 5) + "_lut6.blif");
			fmt::print("mapped network stats    gates {:7d} levels {:7d}\n", mapped.size(), depth(mapped));
			return;
		}

		prev_size = ntk.num_gates();
		prev_levels = depth(ntk);
		mockturtle::write_blif(ntk, outpath);
	}
}

/* runs a node resynthesis iteration on ntk, stops if size does not decrease (for blif)*/
template <typename Ntk, typename Resyn, typename MapFn, typename SaveFn, typename RdFn>
void resyn_iter_blif(Ntk ntk, std::string path, std::string outpath, std::string name, Resyn&& resyn, MapFn&& mapfn, SaveFn&& savefn, RdFn&& read_as_lut, uint32_t maxiter = 1000) {
	uint64_t prev_size = 1000000000u; 
	uint64_t prev_levels = 1000000000u;

	for (auto i = 1u; i <= maxiter; i++) {
		klut_network klut;
		if (i == 1u) {
			fmt::print("read blif as a lut4\n");
			klut = read_as_lut(path, name, 4);
			print_lutmap(klut, i);
		} else {
			klut = mapfn(ntk, path, name, 4);
			print_lutmap(klut, i);
		}

		ntk = node_resynthesis<Ntk>(klut, resyn);
		ntk = cleanup_dangling(ntk);
		print_after(ntk, i);
		
		if (ntk.num_gates() >= prev_size) {
			fmt::print("optimized network stats gates {:7d} levels {:7d}\n", prev_size, prev_levels);
			klut_network mapped = savefn(outpath, outpath.substr(0, outpath.length() - 5) + "_lut6.blif");
			fmt::print("mapped network stats    gates {:7d} levels {:7d}\n", mapped.size(), depth(mapped));
			return;
		}

		prev_size = ntk.num_gates();
		prev_levels = depth(ntk);
		mockturtle::write_blif(ntk, outpath);
	}	
}

#define is_blif(path) (path.length() > 5 && path.substr(path.length() - 5, 5) == ".blif")
#define is_verilog(path) (path.length() > 2 && path.substr(path.length() - 2, 2) == ".v")

/* cut-rewriting with exact mig-npn (built-in) */
void mig_npn_cut_rewrite(std::string path, std::string outpath, int iterations = 1000) {
	fmt::print("mig-npn cut-rewriting [{0}]\n", path);

	mig_network mig;
	mig_npn_resynthesis resyn;

	cut_rewrite_iter_verilog(mig, path, outpath, resyn, iterations); 
}

/* node-resynthesis with exact mig-npn (built-in) */
template <typename MapFn, typename SaveFn, typename RdFn>
void mig_npn_node_resyn(std::string path, std::string outpath, MapFn&& mapfn, SaveFn&& savefn, RdFn&& read_as_lut, int iterations = 1000) {
	fmt::print("mig-npn node-resynthesis [{0}]\n", path);

	mig_npn_resynthesis resyn;
	mig_network mig;

	if (is_verilog(path))
		resyn_iter_verilog(mig, path, outpath, "mig_npn_nr", resyn, mapfn, savefn, iterations);
	else 
		resyn_iter_blif(mig, path, outpath, "mig_npn_nr", resyn, mapfn, savefn, read_as_lut, iterations);
}

/* cut-rewriting with mig (using tig) */
void mig_cut_rewrite(std::string path, std::string outpath, int iterations = 1000) {
	fmt::print("mig cut-rewriting [{0}]\n", path);

	using mig_network = mockturtle::tig_network<three_input_function::majority>;
	mig_network mig;
	
	auto npn_chain = make_shared<tt_chain_map_t>();
	read_npn_to_chain_map(npn_chain, "maj_npn.db");
	exact_tig_resynthesis_params ps{ .cache = npn_chain };
	exact_mig_resynthesis<mig_network> resyn(ps);

	cut_rewrite_iter_verilog(mig, path, outpath, resyn, iterations); 
}


/* node-resynthesis with mig (using tig) */
template <typename MapFn, typename SaveFn, typename RdFn>
void mig_node_resyn(std::string path, std::string outpath, MapFn&& mapfn, SaveFn&& savefn, RdFn&& read_as_lut, int iterations = 1000) {
	fmt::print("mig node-resynthesis [{0}]\n", path);

	using mig_network = mockturtle::tig_network<three_input_function::majority>;
	mig_network mig;

	auto npn_chain = make_shared<tt_chain_map_t>();
	read_npn_to_chain_map(npn_chain, "maj_npn.db");
	exact_tig_resynthesis_params ps{ .cache = npn_chain };
	exact_mig_resynthesis<mig_network> resyn(ps);

	if (is_verilog(path))
		resyn_iter_verilog(mig, path, outpath, "mig_nr", resyn, mapfn, savefn, iterations);
	else
		resyn_iter_blif(mig, path, outpath, "mig_nr", resyn, mapfn, savefn, read_as_lut, iterations);
}

/* cut-rewriting with dig */
void dig_cut_rewrite(std::string path, std::string outpath, int iterations = 1000) {
	fmt::print("dig cut-rewriting [{0}]\n", path);

	using dig_network = mockturtle::tig_network<three_input_function::dot>;
	dig_network dig;

	auto npn_chain = make_shared<tt_chain_map_t>();
	read_npn_to_chain_map(npn_chain, "dot_npn.db");
	exact_tig_resynthesis_params ps{ .cache = npn_chain };
	exact_dig_resynthesis<dig_network> resyn(ps);

	cut_rewrite_iter_verilog(dig, path, outpath, resyn, iterations);
}



/* node-resynthesis with dig */
template <typename MapFn, typename SaveFn, typename RdFn>
void dig_node_resyn(std::string path, std::string outpath,  MapFn&& mapfn, SaveFn&& savefn, RdFn&& read_as_lut, int iterations = 1000) {
	fmt::print("dig node-resynthesis [{0}]\n", path);

	using dig_network = mockturtle::tig_network<three_input_function::dot>;
	dig_network dig;

	auto npn_chain = make_shared<tt_chain_map_t>();
	read_npn_to_chain_map(npn_chain, "dot_npn.db");
	exact_tig_resynthesis_params ps{ .cache = npn_chain };
	exact_dig_resynthesis<dig_network> resyn(ps);

	if (is_verilog(path))
	resyn_iter_verilog(dig, path, outpath, "dig_nr", resyn, mapfn, savefn, iterations);
	else
	resyn_iter_blif(dig, path, outpath, "dig_nr", resyn, mapfn, savefn, read_as_lut, iterations);
}

/* cut-rewriting with aig */
void aig_cut_rewrite(std::string path, std::string outpath, int iterations = 1000) {
	fmt::print("aig cut-rewriting [{0}]\n", path);

	aig_network aig;

	auto npn_cache = make_shared<exact_resynthesis_params::cache_map_t>();
	read_npn_to_chain_map_2(npn_cache, "and_npn.db");
	exact_resynthesis_params ps{ .cache = npn_cache, .blacklist_cache = nullptr};
	exact_aig_resynthesis resyn(false, ps);

	cut_rewrite_iter_verilog(aig, path, outpath, resyn, iterations);
}

/* node-resynthesis with aig */
template <typename MapFn, typename SaveFn, typename RdFn>
void aig_node_resyn(std::string path, std::string outpath, MapFn&& mapfn, SaveFn&& savefn, RdFn&& read_as_lut, int iterations = 1000) {
	fmt::print("aig node-resynthesis [{0}]\n", path);

	aig_network aig;

	auto npn_cache = make_shared<exact_resynthesis_params::cache_map_t>();
	read_npn_to_chain_map_2(npn_cache, "and_npn.db");
	exact_resynthesis_params ps{ .cache = npn_cache, .blacklist_cache = nullptr};
	exact_aig_resynthesis resyn(false, ps);

	if (is_verilog(path))
	resyn_iter_verilog(aig, path, outpath, "aig_nr", resyn, mapfn, savefn, iterations);
	else
	resyn_iter_blif(aig, path, outpath, "aig_nr", resyn, mapfn, savefn, read_as_lut, iterations);
}


template <three_input_function NodeFunc>
void tig_cut_rewrite(std::string path, std::string outpath, std::string dbpath, int iterations = 1000) {
	fmt::print("tig cut-rewriting [{0}]\n", path);

	using tig_network = mockturtle::tig_network<NodeFunc>;
	tig_network tig;

	auto npn_chain = make_shared<tt_chain_map_t>();
	read_npn_to_chain_map(npn_chain, dbpath);
	exact_tig_resynthesis_params ps{ .cache = npn_chain };
	exact_tig_resynthesis<NodeFunc, tig_network> resyn(ps);

	cut_rewrite_iter_verilog(tig, path, outpath, resyn, iterations);
}

template <three_input_function NodeFunc, typename MapFn, typename SaveFn, typename RdFn>
void tig_node_resyn(std::string path, std::string outpath,  std::string dbpath, std::string name, MapFn&& mapfn, SaveFn&& savefn, RdFn&& read_as_lut, int iterations = 1000) {
	fmt::print("tig node-resynthesis [{0}]\n", path);

	using tig_network = mockturtle::tig_network<NodeFunc>;
	tig_network tig;

	auto npn_chain = make_shared<tt_chain_map_t>();
	read_npn_to_chain_map(npn_chain, dbpath);
	exact_tig_resynthesis_params ps{ .cache = npn_chain };
	exact_tig_resynthesis<NodeFunc, tig_network> resyn(ps);

	if (is_verilog(path))
		resyn_iter_verilog(tig, path, outpath, name, resyn, mapfn, savefn, iterations);
	else
		resyn_iter_blif(tig, path, outpath, name, resyn, mapfn, savefn, read_as_lut, iterations);
}

int main(int argc, char ** argv)
{
	(void) argc;

	if (argc <= 1) {
		std::string exec(argv[0]);
		fmt::print("Usage: {} command [options]\n", exec); 
	}

	std::string cmd(argv[1]);

	/* synthesize an AIG for a given function and write the chain to stdout */
	if (cmd == "exact_synth_and2") { 
		if (argc > 2) {
			std::string hex(argv[2]);
			kitty::dynamic_truth_table dt(4);
			kitty::create_from_hex_string(dt, hex);
			synthesize_func_and2(dt);
		} else {
			fmt::print("Truthtable not provided!\n");
		}
		return 0;
	}

	/* generate the npn-to-chain database for a given function */
	std::map<std::string,std::function<void(void)>> m = {
		{"dot", create_dot_database},
		{"onehot", create_onehot_database},
		{"mux", create_mux_database},
		{"andxor", create_andxor_database},
		{"xorand", create_xorand_database},
		{"gamble", create_gamble_database},
		{"orand", create_orand_database},
		{"majority", create_maj_database},
		{"and3", create_and3_database}
	};
	
	if (cmd == "gen_db") {
		if (argc >= 3) {
			std::string gate(argv[2]);
			if (m.count(gate)) {
				m[gate]();
			} else {
				fmt::print("Unsupported gate type!\n");
			}
	
		} else {
			fmt::print("Gate type not provided!\n");
		}
		return 0;
	}

	/* Run experiments */
	std::string exp(argv[1]);
	std::string path(argv[2]);
	std::string outpath(argv[3]);

	int iter = 1000;
	if (argc > 4) {
		std::string iter_str(argv[4]);
		iter = std::stoi(iter_str);
	}

//	if (exp == "mig_npn_cr") {
//		if (is_verilog(path)) {
//			mig_npn_cut_rewrite(path, outpath, iter);
//		} else {
//			std::cerr << "Expected a verilog file for cut-rewriting\n";
//		}
//		return 0;
//	} 
//	if (exp == "mig_npn_nr_old") {
//		mig_npn_node_resyn(path, outpath, lut_map_old<mig_network>, save_mapped_old, read_as_lut4_old, iter);
//		return 0;
//	} 
//	if (exp == "mig_npn_nr_new") {
//		mig_npn_node_resyn(path, outpath, lut_map_new<mig_network>, save_mapped_new, read_as_lut4_new, iter);
//		return 0;
//	} 

	if (exp == "mig_nr") {
		tig_node_resyn<three_input_function::majority>(path, outpath, "maj_npn.db", "mig_nr", 
				lut_map_old<tig_network<three_input_function::majority>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	} 

	if (exp == "dig_nr") {
		tig_node_resyn<three_input_function::dot>(path, outpath, "dot_npn.db", "dig_nr", 
				lut_map_old<tig_network<three_input_function::dot>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	} 	
	
	if (exp == "ohig_nr") {
		tig_node_resyn<three_input_function::onehot>(path, outpath, "onehot_npn.db", "ohig_nr", 
				lut_map_old<tig_network<three_input_function::onehot>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	}

	if (exp == "axig_nr") {
		tig_node_resyn<three_input_function::andxor>(path, outpath, "andxor_npn.db", "axig_nr", 
				lut_map_old<tig_network<three_input_function::andxor>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	}

	if (exp == "xaig_nr") {
		tig_node_resyn<three_input_function::xorand>(path, outpath, "xorand_npn.db", "xaig_nr", 
				lut_map_old<tig_network<three_input_function::xorand>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	}

	if (exp == "muxig_nr") {
		tig_node_resyn<three_input_function::mux>(path, outpath, "mux_npn.db", "muxig_nr", 
				lut_map_old<tig_network<three_input_function::mux>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	}

	if (exp == "oaig_nr") {
		tig_node_resyn<three_input_function::orand>(path, outpath, "orand_npn.db", "oaig_nr", 
				lut_map_old<tig_network<three_input_function::orand>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	}

	if (exp == "gig_nr") {
		tig_node_resyn<three_input_function::gamble>(path, outpath, "gamble_npn.db", "gig_nr", 
				lut_map_old<tig_network<three_input_function::gamble>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	} 


	
	if (exp == "mig_cr") {
		if (is_verilog(path)) {
			tig_cut_rewrite<three_input_function::majority>(path, "new_" + outpath, "maj_npn.db", iter);
		} else {
			std::cerr << "Expected a verilog file for cut-rewriting\n";
		}
		return 0;
	} 








	if (exp == "mig_nr_old") {
		mig_node_resyn(path, outpath, lut_map_old<tig_network<three_input_function::majority>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	} 
	if (exp == "mig_nr_new") {
		mig_node_resyn(path, outpath, lut_map_new<tig_network<three_input_function::majority>>, save_mapped_new, read_as_lut4_new, iter);
		return 0;
	} 



	if (exp == "dig_cr") {
		if (is_verilog(path)) {
			dig_cut_rewrite(path, outpath, iter);
			tig_cut_rewrite<three_input_function::dot>(path, "new_" + outpath, "dot_npn.db", iter);
		} else {
			std::cerr << "Expected a verilog file for cut-rewriting\n";
		}
		return 0;
	}
	if (exp == "ohig_cr") {
		if (is_verilog(path)) {
			tig_cut_rewrite<three_input_function::onehot>(path, "new_" + outpath, "onehot_npn.db", iter);
		} else {
			std::cerr << "Expected a verilog file for cut-rewriting\n";
		}
		return 0;
	}
	if (exp == "gig_cr") {
		if (is_verilog(path)) {
			tig_cut_rewrite<three_input_function::gamble>(path, "new_" + outpath, "gamble_npn.db", iter);
		} else {
			std::cerr << "Expected a verilog file for cut-rewriting\n";
		}
		return 0;
	}
	if (exp == "dig_nr_old") {
		dig_node_resyn(path, outpath, lut_map_old<tig_network<three_input_function::dot>>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	} 
	if (exp == "dig_nr_new") {
		dig_node_resyn(path, outpath, lut_map_new<tig_network<three_input_function::dot>>, save_mapped_new, read_as_lut4_new, iter);
		return 0;
	}  
	
	

	if (exp == "aig_cr") {
		if (is_verilog(path)) {
			aig_cut_rewrite(path, outpath, iter);
		} else {
			std::cerr << "Expected a verilog file for cut-rewriting\n";
		}
		return 0;
	} 
	if (exp == "aig_nr_old") {
		aig_node_resyn(path, outpath, lut_map_old<aig_network>, save_mapped_old, read_as_lut4_old, iter);
		return 0;
	} 
	if (exp == "aig_nr_new") {
		aig_node_resyn(path, outpath, lut_map_new<aig_network>, save_mapped_new, read_as_lut4_new, iter);
		return 0;
	} 
	
	

	std::cerr << "Unrecognized command\n";
	return 0;
}
