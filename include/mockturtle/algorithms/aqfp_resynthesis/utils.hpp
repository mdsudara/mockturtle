#pragma once

#include <functional>
#include <map>
#include <set>
#include <vector>

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
      for ( auto& y : x )
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
