#include <map>
#include <cmath>
#include <ostream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <simple_cfd/cse/cube.hpp>

namespace cfd
{

namespace cse
{

void Cube::incrementUseCounts(std::map<unsigned, std::size_t>& freqs) const
{
  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *literalExponents)
  {
    ++freqs[lMapping.first];
  }
}

int Cube::minExponent(const int a, const int b)
{
  const bool aSign = a < 0;
  const bool bSign = b < 0;

  if (aSign != bSign)
    return 0;

  return (std::abs(a) < std::abs(b)) ? a : b;
}

Cube& Cube::merge(const Cube& c, const bool negate)
{
  exponent_map_t::iterator expIter = begin();

  BOOST_FOREACH(const exponent_map_t::value_type& clMapping, *c.literalExponents)
  {
    const exponent_map_t::value_type 
      lMapping(clMapping.first, negate ? -clMapping.second : clMapping.second);

    while(expIter != literalExponents->end() && expIter->first < lMapping.first)
      ++expIter;

    if (expIter == literalExponents->end() || expIter->first != lMapping.first)
    {
      literalExponents->insert(expIter, lMapping);
    }
    else
    {
      expIter->second += lMapping.second;

      if (expIter->second == 0)
      {
        const exponent_map_t::iterator nextExpIter(boost::next(expIter));
        literalExponents->erase(expIter);
        expIter = nextExpIter;
      }
    }
  }
  return *this;
}

bool Cube::contains(const Cube& c) const
{
  exponent_map_t::const_iterator expIter = literalExponents->begin();

  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *c.literalExponents)
  {
    while(expIter != literalExponents->end() && expIter->first < lMapping.first)
      ++expIter;

    // If c contains a literal not present in the current Cube
    if (expIter == literalExponents->end() 
        || expIter->first != lMapping.first)
      return false;

    const int exponent = minExponent(expIter->second, lMapping.second);

    if (exponent != lMapping.second)
      return false;
  }
  return true;
}

Cube Cube::operator-() const
{
  Cube result(*this);
  BOOST_FOREACH(exponent_map_t::value_type& lMapping, *result.literalExponents)
  {
    lMapping.second = -lMapping.second;
  }
  return result;
}

Cube& Cube::operator+=(const Cube& c)
{
  return merge(c, false);
}

Cube& Cube::operator-=(const Cube& c)
{
  return merge(c, true);
}

Cube& Cube::operator&=(const Cube& c)
{
  exponent_map_t::iterator expIter = literalExponents->begin();

  while(expIter != literalExponents->end())
  {
    const exponent_map_t::iterator nextExpIter(boost::next(expIter));
    const exponent_map_t::const_iterator cIter(c.literalExponents->find(expIter->first));

    if (cIter == c.literalExponents->end()
        || minExponent(expIter->second, cIter->second) == 0)
    {
      literalExponents->erase(expIter);
    }
    else
    {
      expIter->second = minExponent(expIter->second, cIter->second);
    }

    expIter = nextExpIter;
  }
  return *this;
}

std::size_t Cube::numMultiplies() const
{
  std::size_t result = 0;
  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *literalExponents)
  {
    result += std::abs(lMapping.second);
  }
  return result > 0 ? result-1 : 0;
}

std::ostream& operator<<(std::ostream& o, const Cube& c)
{
  c.write(o, detail::LiteralWriter());
  return o;
}

Cube merge(const Cube& a, const Cube& b, const Cube& c)
{
  Cube result(a);
  result+=b;
  result+=c;

  return result;
}

}

}
