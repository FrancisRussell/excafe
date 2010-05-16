#ifndef SIMPLE_CFD_NUMERIC_MONOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_MONOMIAL_HPP

#include <map>
#include <cstddef>
#include <cmath>
#include <set>
#include <iosfwd>
#include <utility>
#include <ostream>
#include <boost/foreach.hpp>
#include <boost/operators.hpp>
#include <simple_cfd/util/lazy_copy.hpp>
#include "value_map.hpp"
#include "cast.hpp"

namespace cfd
{

template<typename V, typename T>
class Monomial : boost::multipliable<Monomial<V,T>,
                 boost::totally_ordered<Monomial<V,T>
                 > >
{
public:
  typedef V variable_t;
  typedef T value_type;

private:
  typedef std::map<variable_t, std::size_t> exponent_map_t;
  util::LazyCopy<exponent_map_t> exponents;

public:
  Monomial()
  {
  }

  Monomial(const Monomial& m) : exponents(m.exponents)
  {
  }

  Monomial(const variable_t& variable, const std::size_t exponent)
  {
    exponents->insert(std::make_pair(variable, exponent));
  }

  Monomial& operator*=(const Monomial& m)
  {
    typename exponent_map_t::iterator expIter = exponents->begin();

    BOOST_FOREACH(const typename exponent_map_t::value_type& mMapping, *m.exponents)
    {
      while(expIter != exponents->end() && expIter->first < mMapping.first)
        ++expIter;

      if (expIter == exponents->end() || !(expIter->first == mMapping.first))
        exponents->insert(expIter, mMapping);
      else
        expIter->second += mMapping.second;
    }

    return *this;
  }

  bool operator==(const Monomial& m) const
  {
    return *exponents == *m.exponents;
  }

  bool operator<(const Monomial& m) const
  {
    return *exponents < *m.exponents;
  }

  Monomial& operator=(const Monomial& m)
  {
    exponents = m.exponents;
    return *this;
  }

  bool isOne() const
  {
    return exponents->empty();
  }

  std::set<variable_t> getVariables() const
  {
    std::set<variable_t> variables;

    for (typename std::map<variable_t, std::size_t>::const_iterator expIter(exponents->begin()); expIter!=exponents->end(); ++expIter)
      variables.insert(expIter->first);

    return variables;
  }

  std::size_t getExponent(const variable_t& variable) const
  {
    const typename std::map<variable_t, std::size_t>::const_iterator varIter(exponents->find(variable));
  
    if (varIter == exponents->end())
    {
      return 0;
    }
    else
    {
      return varIter->second;
    }
  }

  std::pair<value_type, Monomial> derivative(const variable_t& variable) const
  {
    const value_type coefficient = cfd::numeric_cast<double>(getExponent(variable));
    Monomial result;

    for (typename std::map<variable_t, std::size_t>::const_iterator eIter(exponents->begin()); eIter!=exponents->end(); ++eIter)
    {
      if (eIter->first != variable)
      {
        result.exponents->insert(*eIter);
      }
      else if (eIter->second > 1)
      {
        result.exponents->insert(std::make_pair(eIter->first, eIter->second - 1));
      }
    }

    return std::make_pair(coefficient, result);
  }

  std::pair<value_type, Monomial> substituteValues(const detail::ValueMap<variable_t, value_type>& valueMap) const
  {
    Monomial result(*this);
    value_type coefficient = 1.0;

    typedef std::pair<variable_t, value_type> pair_t;

    BOOST_FOREACH(const pair_t& mapping, valueMap.getReference())
    {
      const typename std::map<variable_t, std::size_t>::iterator varIter(result.exponents->find(mapping.first));

      if (varIter != result.exponents->end())
      {
        coefficient *= pow(mapping.second, varIter->second);
        result.exponents->erase(varIter);
      }
    }

    return std::make_pair(coefficient, result);
  }

  void swap(Monomial& m)
  {
    exponents.swap(m.exponents); 
  }

  std::ostream& write(std::ostream& out) const
  {
    for (typename std::map<variable_t, std::size_t>::const_iterator eIter(exponents->begin()); eIter!=exponents->end(); ++eIter)
    {
      out << eIter->first;

      if (eIter->second != 1)
        out << "^{" << eIter->second << "}";
    }
  
    if (exponents->empty())
      out << "1.0";

    return out;
  }

  Monomial substitute(const variable_t& from, const variable_t& to) const
  {
    Monomial result;

    BOOST_FOREACH(const typename exponent_map_t::value_type& exponentMapping, exponents.cref())
    {
      if (exponentMapping.first == from)
        result.exponents->insert(std::make_pair(to, exponentMapping.second));
      else
        result.exponents->insert(exponentMapping);
    }

    return result;
  }
};

template<typename V, typename T>
std::ostream& operator<<(std::ostream& o, const Monomial<V,T>& p)
{
  return p.write(o);
}

}

namespace std
{
  template<typename V, typename T>
  void swap(cfd::Monomial<V,T>& a, cfd::Monomial<V,T>& b)
  {
    a.swap(b);
  }
}

#endif
