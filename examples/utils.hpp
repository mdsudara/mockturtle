#include <cstdio>
#include <percy/percy.hpp>
#include <chrono>
#include <set>
#include <map>
#include <string>
#include <future>
#include <fmt/format.h>

#define MAX_TESTS 256

using namespace percy;
using namespace std;
using namespace mockturtle;

using kitty::dynamic_truth_table;

using tt_chain_map_t = std::unordered_map<kitty::dynamic_truth_table, percy::chain, kitty::hash<kitty::dynamic_truth_table>>;
using tt_tt_map_t = std::unordered_map<kitty::dynamic_truth_table, kitty::dynamic_truth_table, kitty::hash<kitty::dynamic_truth_table>>;

using tt_chain_map_ptr = std::shared_ptr<tt_chain_map_t>;
using tt_tt_map_ptr = std::shared_ptr<tt_tt_map_t>;

template <typename TT>
inline TT ternary_and(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return a & b & c; });
}

template <typename TT>
inline TT ternary_andxor(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return a ^ (b & c); });
}

template <typename TT>
inline TT ternary_dot(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return a ^ (c | (a & b)); });
}

template <typename TT>
inline TT ternary_gamble(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return ~((a & b & c) ^ (~a & ~b & ~c)); });
}

template <typename TT>
inline TT ternary_majority(const TT &first, const TT &second, const TT &third)
{
    return kitty::ternary_majority(first, second, third);
}

template <typename TT>
inline TT ternary_mux(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return (a & b) | (~a & c); });
}

template <typename TT>
inline TT ternary_onehot(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return (a & ~b & ~c) ^ (~a & b & ~c) ^ (~a & ~b & c); });
}

template <typename TT>
inline TT ternary_orand(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return a & (b | c); });
}

template <typename TT>
inline TT ternary_xor(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return a ^ b ^ c; });
}

template <typename TT>
inline TT ternary_xorand(const TT &first, const TT &second, const TT &third)
{
    return ternary_operation(first, second, third, [](auto a, auto b, auto c) { return a & (b ^ c); });
}

template <typename Fn>
void generate_code_for_primitive_mapping(Fn &&fn)
{
    kitty::dynamic_truth_table a{3};
    kitty::dynamic_truth_table b{3};
    kitty::dynamic_truth_table c{3};
    kitty::dynamic_truth_table cnst0{3};
    
    kitty::create_nth_var(a, 0);
    kitty::create_nth_var(b, 1);
    kitty::create_nth_var(c, 2);
    kitty::create_from_hex_string(cnst0, "00");
    
    std::array<kitty::dynamic_truth_table, 8> choices = {a, ~a, b, ~b, c, ~c, cnst0, ~cnst0};
    
    std::set<std::string> all;

    for (auto ind0 = 0; ind0 < 8; ind0++)
    {
        for (auto ind1 = 0; ind1 < 8; ind1++)
        {
            for (auto ind2 = 0; ind2 < 8; ind2++)
            {
                kitty::dynamic_truth_table prim = fn(choices[ind0], choices[ind1], choices[ind2]);
                std::array<int, 3> inds = {ind0, ind1, ind2};
          		
				std::string code = "case 0x{4}:\n\tsignals.emplace_back( {5}ntk.create_dot( "; 

                for (int j = 0; j < 3; j++)
                {
                    switch (inds[j])
                    {
                        case 0:
							code += "{0}";
                            break;
                        case 1:
                            code += "!{0}";
                            break;
                        case 2:
                            code += "{1}";
                            break;
                        case 3:
                            code += "!{1}";
                            break;
                        case 4:
                            code += "{2}";
                            break;
                        case 5:
                            code += "!{2}";
                            break;
                        case 6:
                            code += "{3}";
                            break;
                        case 7:
                            code += "!{3}";
                            break;
                    }
                    if (j < 2) {
						code += ", ";
					}
						
                }
                code += " ) );\n\tbreak;\n";

				std::string output_inv = "";
                
                if (!kitty::is_normal(prim))
                {
                    prim = ~prim;
					output_inv = "!";
                }
                
                auto ret = all.insert(kitty::to_hex(prim));
                if (ret.second)
                {
					fmt::print(code, "c1", "c2", "c3", "ntk.get_constant( false )", kitty::to_hex(prim), output_inv);
                }
            }
        }
    }
}

template <typename Fn>
void add_all_three_input_primitives(percy::spec &spec, Fn &&fn)
{
    
    kitty::dynamic_truth_table a{3};
    kitty::dynamic_truth_table b{3};
    kitty::dynamic_truth_table c{3};
    kitty::dynamic_truth_table cnst0{3};
    
    kitty::create_nth_var(a, 0);
    kitty::create_nth_var(b, 1);
    kitty::create_nth_var(c, 2);
    kitty::create_from_hex_string(cnst0, "00");
    
    std::array<kitty::dynamic_truth_table, 8> choices = {a, ~a, b, ~b, c, ~c, cnst0, ~cnst0};
    
    std::set<std::string> all;
    
    for (auto ind0 = 0; ind0 < 8; ind0++)
    {
        for (auto ind1 = 0; ind1 < 8; ind1++)
        {
            for (auto ind2 = 0; ind2 < 8; ind2++)
            {
                kitty::dynamic_truth_table prim = fn(choices[ind0], choices[ind1], choices[ind2]);
                
                if (!kitty::is_normal(prim))
                {
                    prim = ~prim;
                }
                
                auto ret = all.insert(kitty::to_hex(prim));
                if (ret.second)
                {
                    spec.add_primitive(prim);
                }
            }
        }
    }
}

template <typename Fn>
percy::chain synthesize_func_general_primitive(kitty::dynamic_truth_table const &tt, Fn &&fn)
{
    percy::chain chain;
    percy::spec spec;

    spec.fanin = 3;
    spec[0] = tt;

    add_all_three_input_primitives(spec, fn);
    
    auto result = percy::synthesize(spec, chain);

    assert(result == percy::success);
    assert(chain.simulate()[0] == spec[0]);
    
    return chain;
}

template <typename Fn>
void synthesize_all_npn_classes(Fn &&fn, int num_inputs, tt_chain_map_ptr npn_chain)
{
	std::set<kitty::dynamic_truth_table> npn;
    
    kitty::dynamic_truth_table tt(num_inputs);
    do
    {
        const auto res = std::get<0>(kitty::exact_npn_canonization(tt));
		npn.insert(res);

        kitty::next_inplace(tt);
    } while (!kitty::is_const0(tt));

    for (auto it = npn.begin(); it != npn.end(); it ++)
    {
		if (kitty::is_normal(*it)) {
       		(*npn_chain)[*it] = synthesize_func_general_primitive(*it, fn);
		} else {
       		(*npn_chain)[~(*it)] = synthesize_func_general_primitive(~(*it), fn);
		}
    }
}


percy::chain synthesize_func_and2(kitty::dynamic_truth_table const &tt)
{
	percy::chain chain;
	percy::spec spec;

	spec.fanin = 2;
	spec[0] = tt;

	spec.set_primitive( percy::AIG );
	spec.verbosity = 0;
	spec.add_alonce_clauses = true;
	spec.add_colex_clauses = true; 
	spec.add_lex_clauses = false; 
	spec.add_lex_func_clauses = true;
	spec.add_nontriv_clauses = true; 
	spec.add_noreapply_clauses = true; 
	spec.add_symvar_clauses = true; 
	spec.conflict_limit = 0; 

	//std::shared_ptr<percy::glucose_wrapper> solver = std::make_shared<percy::glucose_wrapper>();
	//auto encoder = get_encoder(*solver, percy::ENC_SSV);
	//solver->set_nr_threads(96);
	auto result = percy::synthesize(spec, chain); // *solver, *encoder);

	assert(result == percy::success);
	assert(chain.simulate()[0] == spec[0]);

	std::cout << kitty::to_hex(tt) << " " << chain.get_nr_steps() << " " << chain.get_outputs()[0];
	for ( auto i = 0; i < chain.get_nr_steps(); i++ ) {
		std::cout << " " << kitty::to_hex( chain.get_operator(i) );
		std::cout << " " << chain.get_step(i)[0] << " " << chain.get_step(i)[1];
	}
	std::cout << std::endl;

	return chain;
}


void synthesize_all_npn_classes_2(int num_inputs, tt_chain_map_ptr npn_chain)
{
	std::set<kitty::dynamic_truth_table> npn;
    
    kitty::dynamic_truth_table tt(num_inputs);
    do
    {
        const auto res = std::get<0>(kitty::exact_npn_canonization(tt));
		npn.insert(res);

        kitty::next_inplace(tt);
    } while (!kitty::is_const0(tt));

	auto count = 0u;
    for (auto it = npn.begin(); it != npn.end(); it ++)
    {
		fmt::print("Synthesizing {0} out of {1}: TT = {2}\n", ++count, npn.size(), kitty::to_hex(*it));
		if (kitty::is_normal(*it)) {
       		(*npn_chain)[*it] = synthesize_func_and2(*it);
		} else {
       		(*npn_chain)[~(*it)] = synthesize_func_and2(~(*it));
		}
    }
}


void read_npn_to_chain_map( tt_chain_map_ptr cache, std::string path ) {
	std::ifstream in(path);

	std::string tt_hex, op_hex;
	int chain_len, output, i1, i2, i3;

	for (auto i = 0; i < 222; i++) {	
		in >> tt_hex;

		percy::chain ch;
		in >> chain_len;
		in >> output;

		ch.reset(4, 1, chain_len, 3);
		ch.set_output(0, output);

		for (auto i = 0; i < chain_len; i++) {
			in >> op_hex;
			kitty::dynamic_truth_table dop(3);
			kitty::create_from_hex_string(dop, op_hex);
			in >> i1 >> i2 >> i3;
			std::vector inputs{i1, i2, i3};
			ch.set_step(i, inputs, dop);
		}

		kitty::dynamic_truth_table dtt(4);
		kitty::create_from_hex_string(dtt, tt_hex);
		(*cache)[dtt] = ch;
	}

	in.close();
}

void read_npn_to_chain_map_2( tt_chain_map_ptr cache, std::string path ) {
	std::ifstream in(path);

	std::string tt_hex, op_hex;
	int chain_len, output, i1, i2;

	for (auto i = 0; i < 222; i++) {	
		in >> tt_hex;

		percy::chain ch;
		in >> chain_len;
		in >> output;

		ch.reset(4, 1, chain_len, 2);
		ch.set_output(0, output);

		for (auto i = 0; i < chain_len; i++) {
			in >> op_hex;
			kitty::dynamic_truth_table dop(2);
			kitty::create_from_hex_string(dop, op_hex);
			in >> i1 >> i2;
			std::vector inputs{i1, i2};
			ch.set_step(i, inputs, dop);
		}

		kitty::dynamic_truth_table dtt(4);
		kitty::create_from_hex_string(dtt, tt_hex);
		(*cache)[dtt] = ch;
	}

	in.close();
}

void read_tt_to_npn_map( tt_tt_map_ptr cache, std::string path ) {
	std::ifstream in(path);
	std::string tt_hex, npn_hex;

	for (auto i = 0; i < (1 << 16); i++) {
		in >> tt_hex >> npn_hex;
		kitty::dynamic_truth_table tt(4);
		kitty::dynamic_truth_table npn(4);

		kitty::create_from_hex_string(tt, tt_hex);
		kitty::create_from_hex_string(npn, npn_hex);
		(*cache)[tt] = npn;
	}

	in.close();
}


void write_npn_to_chain_map( tt_chain_map_ptr cache, const std::string& path ) {
	std::ofstream out( path );

	for ( auto it = cache->begin(); it != cache->end(); it++ ) {
		out << kitty::to_hex(it->first) << " " << it->second.get_nr_steps() << " " << it->second.get_outputs()[0];
		for ( auto i = 0; i < it->second.get_nr_steps(); i++ ) {
			out << " " << kitty::to_hex( it->second.get_operator(i) );
			out << " " << it->second.get_step(i)[0] << " " << it->second.get_step(i)[1] << " " << it->second.get_step(i)[2];
		}
		out << std::endl;
	}

	out.close();
}

void write_npn_to_chain_map_2( tt_chain_map_ptr cache, const std::string& path ) {
	std::ofstream out( path );

	for ( auto it = cache->begin(); it != cache->end(); it++ ) {
		out << kitty::to_hex(it->first) << " " << it->second.get_nr_steps() << " " << it->second.get_outputs()[0];
		for ( auto i = 0; i < it->second.get_nr_steps(); i++ ) {
			out << " " << kitty::to_hex( it->second.get_operator(i) );
			out << " " << it->second.get_step(i)[0] << " " << it->second.get_step(i)[1];
		}
		out << std::endl;
	}

	out.close();
}


void write_tt_to_npn( tt_tt_map_ptr cache, const std::string& path ) {
	std::ofstream out( path );
	
	for ( auto it = cache->begin(); it != cache->end(); it++ ) {
		out << kitty::to_hex(it->first) << " " << kitty::to_hex(it->second) << std::endl;
	}
	
	out.close();
}

void create_maj_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(kitty::ternary_majority<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "maj_npn.db");
}


void create_dot_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_dot<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "dot_npn.db");
}

void create_onehot_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_onehot<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "onehot_npn.db");
}

void create_gamble_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_gamble<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "gamble_npn.db");
}

void create_mux_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_mux<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "mux_npn.db");
}

void create_and3_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_and<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "and3_npn.db");
}

void create_orand_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_orand<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "orand_npn.db");
}

void create_xorand_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_xorand<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "xorand_npn.db");
}

void create_andxor_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes(ternary_andxor<kitty::dynamic_truth_table>, 4, npn_chain);
	write_npn_to_chain_map(npn_chain, "andxor_npn.db");
}

void create_and_database() {
	tt_chain_map_ptr npn_chain = make_shared<tt_chain_map_t>();
	synthesize_all_npn_classes_2(4, npn_chain);
	write_npn_to_chain_map_2(npn_chain, "and_npn.db");
}

template<typename Ntk>
mockturtle::klut_network lut_map_old( Ntk const& ntk, std::string design, std::string command, uint32_t k = 4 )
{
	std::string tempfile1 = design + ".temp1old." + command + ".blif";
	std::string tempfile2 = design + ".temp2old." + command + ".blif";

	mockturtle::write_blif( ntk, tempfile1 );
	system( fmt::format( "abc -q \"{}; if -a -K {}; write_blif {}\"", tempfile1, k, tempfile2 ).c_str() );

	mockturtle::klut_network klut;
	if ( lorina::read_blif( tempfile2, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
	{
		std::cout << "FATAL OLD LUT MAP - Reading mapped network failed!" << std::endl;
		std::abort();
		return klut;
	}

	system (  fmt::format("rm {}", tempfile1).c_str() );
	system (  fmt::format("rm {}", tempfile2).c_str() );
	return klut;
}

template<typename Ntk>
mockturtle::klut_network lut_map_new( Ntk const& ntk, std::string design, std::string command, uint32_t k = 4 )
{
	std::string tempfile1 = design + ".temp1new." + command + ".blif";
	std::string tempfile2 = design + ".temp2new." + command + ".blif";

	mockturtle::write_blif( ntk, tempfile1 );
	system( fmt::format( "abc -q \"{}; &get; &if -a -K {}; &put; write_blif {}\"", tempfile1, k, tempfile2 ).c_str() );

	mockturtle::klut_network klut;
	if ( lorina::read_blif( tempfile2, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
	{
		std::cout << "FATAL NEW LUT MAP - Reading mapped network failed!" << std::endl;
		std::abort();
		return klut;
	}

	system (  fmt::format("rm {}", tempfile1).c_str() );
	system (  fmt::format("rm {}", tempfile2).c_str() );
	return klut;
}

mockturtle::klut_network save_mapped_old( std::string src, std::string dst )
{
	system( fmt::format( "abc -q \"{}; if -a -K {}; write_blif {}\"", src, 6, dst ).c_str() );
	mockturtle::klut_network klut;
	lorina::read_blif ( dst, mockturtle::blif_reader( klut ) );
	return klut;
}

mockturtle::klut_network save_mapped_new( std::string src, std::string dst )
{
	system( fmt::format( "abc -q \"{}; &get; &if -a -K {}; &put; write_blif {}\"", src, 6, dst ).c_str() );
	mockturtle::klut_network klut;
	lorina::read_blif ( dst, mockturtle::blif_reader( klut ) );
	return klut;
}

template <typename Ntk>
int depth(Ntk&& ntk) {
	depth_view dview{ntk};
	return dview.depth();
}

mockturtle::klut_network read_as_lut4_old( std::string design, std::string command, uint32_t k = 4 )
{
	std::string tempfile = design + ".tempold." + command + ".blif";
	system( fmt::format( "abc -q \"{}; if -a -K {}; write_blif {}\"", design, k, tempfile ).c_str() );

	mockturtle::klut_network klut;
	if ( lorina::read_blif( tempfile, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
	{
		std::cout << "FATAL - Reading failed!" << std::endl;
		std::abort();
		return klut;
	}

	system (  fmt::format("rm {}", tempfile).c_str() );
	return klut;
}

mockturtle::klut_network read_as_lut4_new( std::string design, std::string command, uint32_t k = 4 )
{
	std::string tempfile = design + ".tempnew." + command + ".blif";
	system( fmt::format( "abc -q \"{}; &get; &if -a -K {}; &put; write_blif {}\"", design, k, tempfile ).c_str() );

	mockturtle::klut_network klut;
	if ( lorina::read_blif( tempfile, mockturtle::blif_reader( klut ) ) != lorina::return_code::success )
	{
		std::cout << "FATAL - Reading failed!" << std::endl;
		std::abort();
		return klut;
	}

	system (  fmt::format("rm {}", tempfile).c_str() );
	return klut;
}

