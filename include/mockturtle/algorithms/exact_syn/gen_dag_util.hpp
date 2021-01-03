#pragma once

#include <chrono>
#include <future>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <vector>

#include <fmt/format.h>

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::set<T>& s );

template<typename T>
std::ostream& operator<<( std::ostream& os, const std::multiset<T>& s );

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

namespace std
{
template<typename A>
struct hash<std::vector<A>>;

template<typename A, typename B, typename C>
struct hash<std::tuple<A, B, C>>;

template<typename A>
struct hash<std::multiset<std::multiset<A>>>;

template<typename A, typename B>
struct hash<std::map<A, B>>;

template<typename A>
struct hash<std::vector<A>>
{
public:
  size_t operator()( const std::vector<A>& key ) const
  {
	std::size_t seed = key.size();
	for ( auto& i : key )
	{
	  seed ^= ha( i ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	}
	return seed;
  }

private:
  std::hash<A> ha;
};

template<typename A, typename B, typename C>
struct hash<std::tuple<A, B, C>>
{
public:
  size_t operator()( const std::tuple<const A, const B, const C>& key ) const
  {
	size_t seed = ha( std::get<0>( key ) );
	seed ^= hb( std::get<1>( key ) ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	seed ^= hc( std::get<2>( key ) ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	return seed;
  }

private:
  std::hash<A> ha;
  std::hash<B> hb;
  std::hash<C> hc;
};

template<typename A>
struct hash<std::multiset<std::multiset<A>>>
{
public:
  size_t operator()( const std::multiset<std::multiset<A>>& key ) const
  {
	size_t seed = key.size();
	for ( auto& x : key )
	{
	  seed ^= x.size() + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	  for ( auto& y : key )
	  {
		seed ^= ha( y ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	  }
	}
	return seed;
  }

private:
  std::hash<A> ha;
};

template<typename A, typename B>
struct hash<std::map<A, B>>
{
public:
  size_t operator()( const std::map<A, B>& key ) const
  {
	size_t seed = key.size();
	for ( auto it = key.begin(); it != key.end(); it++ )
	{
	  seed ^= ha( it->first ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	  seed ^= hb( it->second ) + 0x9e3779b9 + ( seed << 6 ) + ( seed >> 2 );
	}
	return seed;
  }

private:
  std::hash<A> ha;
  std::hash<B> hb;
};

} // namespace std

namespace mockturtle
{

namespace detail
{

/*
 * Computes and returns the frequency map for a given collection of elements
 */
template<typename ElemT>
inline auto get_frequencies( const std::vector<ElemT>& elems )
{
  std::map<ElemT, uint32_t> elem_counts;
  std::for_each( elems.begin(), elems.end(), [&elem_counts]( auto e ) { elem_counts[e]++; } );
  //for ( const auto& e : elems )
  //{
  //elem_counts[e]++;
  //}
  return elem_counts;
}

template<typename ElemT>
class partition_generator
{
  using part = std::multiset<ElemT>;
  using partition = std::multiset<part>;
  using partition_set = std::set<partition>;

  using inner_cache_key_t = std::vector<ElemT>;
  using inner_cache_t = std::unordered_map<inner_cache_key_t, partition_set>;

  using outer_cache_key_t = std::tuple<std::vector<uint32_t>, uint32_t, uint32_t>;
  using outer_cache_t = std::unordered_map<outer_cache_key_t, inner_cache_t>;

public:
  /**
   * \brief Computes and returns a vector of partitions for a given list of elements
   * such that no part contains any element 'e' more than 'max_counts[e]' times.
   */
  partition_set operator()(
	  std::vector<int> elems,
	  const std::vector<uint32_t>& max_counts = {},
	  uint32_t max_parts = 0,
	  uint32_t max_part_size = 0 )
  {
	_elems = elems;
	_max_counts = max_counts;
	_max_parts = max_parts;
	_max_part_size = max_part_size;

	const outer_cache_key_t key = { _max_counts, _max_parts, _max_part_size };
	partition_cache = outer_cache.insert( { key, inner_cache_t() } ).first;

	return get_all_partitions();
  }

private:
  outer_cache_t outer_cache;

  typename outer_cache_t::iterator partition_cache;
  std::vector<ElemT> _elems;
  std::vector<uint32_t> _max_counts;
  uint32_t _max_parts;
  uint32_t _max_part_size;

  partition_set get_all_partitions()
  {
	if ( _elems.size() == 0 )
	{
	  return { {} }; // return the empty partition.
	}

	inner_cache_key_t key = _elems;
	if ( partition_cache->second.count( key ) )
	{
	  return partition_cache->second.at( key );
	}

	partition_set result;

	auto last = _elems.back();
	_elems.pop_back();

	auto temp = get_all_partitions();

	for ( auto&& t : temp )
	{
	  partition cpy;

	  // take 'last' in its own partition
	  cpy = t;

	  if ( _max_parts == 0u || _max_parts > cpy.size() )
	  {
		cpy.insert( { last } );
		result.insert( cpy );
	  }

	  // add 'last' to one of the existing partitions
	  for ( auto it = t.begin(); it != t.end(); )
	  {
		if ( _max_counts.empty() || it->count( last ) < _max_counts[last] )
		{

		  if ( _max_part_size == 0 || _max_part_size > it->size() )
		  {
			cpy = t;
			auto elem_it = cpy.find( *it );
			auto cpy_elem = *elem_it;
			cpy_elem.insert( last );
			cpy.erase( elem_it );
			cpy.insert( cpy_elem );
			result.insert( cpy );
		  }
		}

		std::advance( it, t.count( *it ) );
	  }
	}

	return ( partition_cache->second[key] = result );
  }
};

template<typename ElemT>
class partition_extender
{
  using part = std::multiset<ElemT>;
  using partition = std::multiset<part>;
  using partition_set = std::set<partition>;

  using inner_cache_key_t = std::vector<ElemT>;
  using inner_cache_t = std::unordered_map<inner_cache_key_t, partition_set>;

  using outer_cache_key_t = std::tuple<partition, std::vector<uint32_t>, uint32_t>;
  using outer_cache_t = std::map<outer_cache_key_t, inner_cache_t>;

public:
  /**
   * \brief Compute a list of different partitions that can be obtained by adding elements in
   * 'elems' to the parts of 'base' such that no part contains any element 'e' more than
   * 'max_counts[e]' times
   */
  partition_set operator()( std::vector<ElemT> elems, partition base, const std::vector<uint32_t>& max_counts, uint32_t max_part_size = 0 )
  {
	_elems = elems;
	_base = base;
	_max_counts = max_counts;
	_max_part_size = max_part_size;

	const outer_cache_key_t key = { _base, _max_counts, _max_part_size };
	partition_cache = outer_cache.insert( { key, inner_cache_t() } ).first;

	return extend_partitions();
  }

private:
  outer_cache_t outer_cache;

  typename outer_cache_t::iterator partition_cache;
  std::vector<ElemT> _elems;
  partition _base;
  std::vector<uint32_t> _max_counts;
  uint32_t _max_part_size;

  partition_set extend_partitions()
  {
	if ( _elems.size() == 0 )
	{
	  return { _base };
	}

	inner_cache_key_t key = _elems;
	if ( partition_cache->second.count( key ) )
	{
	  return partition_cache->second.at( key );
	}

	partition_set result;

	auto last = _elems.back();
	_elems.pop_back();

	auto temp = extend_partitions();
	for ( auto&& t : temp )
	{
	  partition cpy;

	  for ( auto it = t.begin(); it != t.end(); )
	  {
		if ( it->count( last ) < _max_counts.at( last ) )
		{

		  if ( _max_part_size == 0 || _max_part_size > it->size() )
		  {
			cpy = t;
			auto elem_it = cpy.find( *it );
			auto cpy_elem = *elem_it;
			cpy_elem.insert( last );
			cpy.erase( elem_it );
			cpy.insert( cpy_elem );
			result.insert( cpy );
		  }
		}

		std::advance( it, t.count( *it ) );
	  }
	}

	return ( partition_cache->second[key] = result );
  }
};

template<typename ElemT>
struct sublist_generator
{
  using sub_list_cache_key_t = std::map<ElemT, uint32_t>;

public:
  /**
   * \brief Given a list of elements 'elems', generate all sub lists of those elements.
   * Ex: if 'elems' = [1, 2, 2, 3], this will generate the following lists:
   * [0], [1], [1, 2], [1, 2, 2], [1, 2, 2, 3], [1, 2, 3], [1, 3], [2], [2, 2], [2, 2, 3], [2, 3], and [3].
   */
  std::set<std::vector<ElemT>> operator()( std::vector<ElemT> elems )
  {
	elem_counts = get_frequencies( elems );
	return get_sub_lists_recur();
  }

private:
  std::unordered_map<sub_list_cache_key_t, std::set<std::vector<ElemT>>> sub_list_cache;
  std::map<ElemT, uint32_t> elem_counts;

  std::set<std::vector<ElemT>> get_sub_lists_recur()
  {
	if ( elem_counts.size() == 0u )
	{
	  return { {} };
	}

	sub_list_cache_key_t key = elem_counts;
	if ( !sub_list_cache.count( key ) )
	{
	  auto last = std::prev( elem_counts.end() );
	  auto last_elem = last->first;
	  auto last_count = last->second;
	  elem_counts.erase( last );

	  std::set<std::vector<int>> result;

	  std::vector<int> t;
	  for ( auto i = last_count; i > 0; --i )
	  {
		t.push_back( last_elem );
		result.insert( t ); // insert a copy of t, and note that t is already sorted.
	  }

	  auto temp = get_sub_lists_recur();

	  for ( std::vector<int> t : temp )
	  {
		result.insert( t );
		for ( auto i = last_count; i > 0; --i )
		{
		  t.push_back( last_elem );
		  std::sort( t.begin(), t.end() );
		  result.insert( t );
		}
	  }

	  sub_list_cache[key] = result;
	}

	return sub_list_cache[key];
  }
};

} // namespace detail

} // namespace mockturtle

