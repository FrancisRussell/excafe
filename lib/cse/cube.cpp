#include <map>
#include <cmath>
#include <ostream>
#include <algorithm>
#include <boost/foreach.hpp>
#include <boost/utility.hpp>
#include <excafe/cse/cube.hpp>
#include <excafe/cse/literal_info.hpp>
#include <excafe/cse/new_literal_creator.hpp>

namespace excafe
{

namespace cse
{

void Cube::incrementUseCounts(std::map<LiteralInfo, std::size_t>& freqs) const
{
  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *literalExponents)
  {
    ++freqs[LiteralInfo(lMapping.first, lMapping.second < 0)];
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

Cube& Cube::merge(const Cube& c2, const bool negate)
{
  typedef exponent_map_t::value_type::second_type exponent_t;
  const Cube& c1 = *this;
  exponent_map_t::const_iterator c1Iter = c1.begin();
  exponent_map_t::const_iterator c2Iter = c2.begin();

  exponent_map_t newExponents;
  newExponents.reserve(c1.size() + c2.size());

  while(c1Iter != c1.end() || c2Iter != c2.end())
  {
    if (c2Iter == c2.end() || (c1Iter != c1.end() && c1Iter->first < c2Iter->first))
    {
      newExponents.insert(*c1Iter);
      ++c1Iter;
    }
    else
    {
      const exponent_t c2NegatedValue = (negate ? -c2Iter->second : c2Iter->second);

      if (c1Iter == c1.end() || c2Iter->first < c1Iter->first)
      {
        newExponents.insert(exponent_map_t::value_type(c2Iter->first, c2NegatedValue));
      }
      else
      {
        const exponent_map_t::value_type::first_type exponentSum = c1Iter->second + c2NegatedValue;
        if (exponentSum != 0)
          newExponents.insert(exponent_map_t::value_type(c1Iter->first, exponentSum));
        ++c1Iter;
      }

      ++c2Iter;
    }
  }

  swap(*this->literalExponents, newExponents);
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
  result *= -1;
  return result;
}

Cube& Cube::operator*=(const int exponent)
{
  BOOST_FOREACH(exponent_map_t::value_type& lMapping, *literalExponents)
  {
    lMapping.second *= exponent;
  }
  return *this;
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
    const exponent_map_t::const_iterator cIter(c.literalExponents->find(expIter->first));

    if (cIter == c.literalExponents->end()
        || minExponent(expIter->second, cIter->second) == 0)
    {
      expIter = literalExponents->erase(expIter);
    }
    else
    {
      expIter->second = minExponent(expIter->second, cIter->second);
      ++expIter;
    }
  }
  return *this;
}

bool Cube::isUnit(const NewLiteralCreator& creator) const
{
  if (literalExponents->empty())
  {
    return true;
  }
  else if (literalExponents->size() == 1 
           && creator.isUnit(literalExponents->begin()->first))
  {
    return true;
  }
  else
  {
    return false;
  }
}

bool Cube::isNumeric(const NewLiteralCreator& creator) const
{
  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *literalExponents)
  {
    if (!creator.isNumeric(lMapping.first))
      return false;
  }

  return true;
}

bool Cube::hasCoefficient(const NewLiteralCreator& creator) const
{
  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *literalExponents)
  {
    if (creator.isNumeric(lMapping.first) && !creator.isUnit(lMapping.first))
      return true;
  }

  return false;
}

std::size_t Cube::numMultiplies(const NewLiteralCreator& creator) const
{
  std::size_t result = 0;
  bool hasNumeric = false;

  BOOST_FOREACH(const exponent_map_t::value_type& lMapping, *literalExponents)
  {
    if (creator.isNumeric(lMapping.first))
      hasNumeric = true;
    else
      result += std::abs(lMapping.second);
  }

  if (hasNumeric)
    ++result;

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
