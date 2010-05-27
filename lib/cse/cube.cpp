#include <map>
#include <ostream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <simple_cfd/util/maybe.hpp>
#include <simple_cfd/cse/cube.hpp>

namespace cfd
{

namespace cse
{

util::Maybe<Cube> Cube::operator/(const Cube& c) const
{
  Cube result(*this);
  std::map<unsigned, unsigned>::iterator expIter = result.literalExponents.begin();

  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, c.literalExponents)
  {
    while(expIter != result.literalExponents.end() && expIter->first < lMapping.first)
      ++expIter;

    if (expIter == result.literalExponents.end() || 
        expIter->first != lMapping.first ||
        lMapping.second > expIter->second)
    {
      return util::Nothing();
    }
    else if (lMapping.second == expIter->second)
    {
      const std::map<unsigned, unsigned>::iterator nextExpIter(boost::next(expIter));
      result.literalExponents.erase(expIter);
      expIter = nextExpIter;
    }
    else if (lMapping.second < expIter->second)
    {
      expIter->second -= lMapping.second;
    }
  }
  return result;
}

Cube& Cube::operator+=(const Cube& c)
{
  std::map<unsigned, unsigned>::iterator expIter = literalExponents.begin();

  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, c.literalExponents)
  {
    while(expIter != literalExponents.end() && expIter->first < lMapping.first)
      ++expIter;

    if (expIter == literalExponents.end() || expIter->first != lMapping.first)
      literalExponents.insert(expIter, lMapping);
    else
      expIter->second += lMapping.second;
  }
  return *this;
}

Cube& Cube::operator&=(const Cube& c)
{
  std::map<unsigned, unsigned>::iterator expIter = literalExponents.begin();

  while(expIter != literalExponents.end())
  {
    const std::map<unsigned, unsigned>::iterator nextExpIter(boost::next(expIter));
    const std::map<unsigned, unsigned>::const_iterator cIter(c.literalExponents.find(expIter->first));

    if (cIter == c.literalExponents.end())
    {
      literalExponents.erase(expIter);
    }
    else
    {
      expIter->second = std::min(expIter->second, cIter->second);
    }

    expIter = nextExpIter;
  }
  return *this;
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
