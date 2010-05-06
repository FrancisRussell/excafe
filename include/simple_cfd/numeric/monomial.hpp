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

namespace cfd
{

template<typename V>
class Monomial : boost::multipliable<Monomial<V> >
{
public:
  typedef V variable_t;

private:
  typedef std::map<variable_t, std::size_t> exponent_map_t;
  exponent_map_t exponents;

public:
  Monomial()
  {
  }

  Monomial(const Monomial& m) : exponents(m.exponents)
  {
  }

  Monomial(const variable_t& variable, const std::size_t exponent)
  {
    exponents.insert(std::make_pair(variable, exponent));
  }

  Monomial& operator*=(const Monomial& m)
  {
    for (typename std::map<variable_t, std::size_t>::const_iterator expIter(m.exponents.begin()); expIter!=m.exponents.end(); ++expIter)
      exponents[expIter->first] += expIter->second;

    return *this;
  }

  bool operator==(const Monomial& m) const
  {
    return exponents == m.exponents;
  }

  bool operator<(const Monomial& m) const
  {
    return exponents < m.exponents;
  }

  Monomial& operator=(const Monomial& m)
  {
    exponents = m.exponents;
    return *this;
  }

  bool isOne() const
  {
    return exponents.empty();
  }

  std::set<variable_t> getVariables() const
  {
    std::set<variable_t> variables;

    for (typename std::map<variable_t, std::size_t>::const_iterator expIter(exponents.begin()); expIter!=exponents.end(); ++expIter)
      variables.insert(expIter->first);

    return variables;
  }

  std::size_t getExponent(const variable_t& variable) const
  {
    const typename std::map<variable_t, std::size_t>::const_iterator varIter(exponents.find(variable));
  
    if (varIter == exponents.end())
    {
      return 0;
    }
    else
    {
      return varIter->second;
    }
  }

  std::pair<double, Monomial> derivative(const variable_t& variable) const
  {
    const double coefficient = getExponent(variable);
    Monomial result;

    for (typename std::map<variable_t, std::size_t>::const_iterator eIter(exponents.begin()); eIter!=exponents.end(); ++eIter)
    {
      if (eIter->first != variable)
      {
        result.exponents.insert(*eIter);
      }
      else if (eIter->second > 1)
      {
        result.exponents.insert(std::make_pair(eIter->first, eIter->second - 1));
      }
    }

    return std::make_pair(coefficient, result);
  }

  std::pair<double, Monomial> substituteValue(const variable_t& var, const double value) const
  {
    Monomial result(*this);
    double coefficient = 1.0;
    const typename std::map<variable_t, std::size_t>::iterator varIter(result.exponents.find(var));

    if (varIter != result.exponents.end())
    {
      coefficient = std::pow(value, varIter->second);
      result.exponents.erase(varIter);
    }

    return std::make_pair(coefficient, result);
  }

  void swap(Monomial& m)
  {
    exponents.swap(m.exponents); 
  }

  std::ostream& write(std::ostream& out) const
  {
    for (typename std::map<variable_t, std::size_t>::const_iterator eIter(exponents.begin()); eIter!=exponents.end(); ++eIter)
    {
      out << eIter->first;

      if (eIter->second != 1)
        out << "^{" << eIter->second << "}";
    }
  
    if (exponents.empty())
      out << "1.0";

    return out;
  }

  Monomial substitute(const variable_t& from, const variable_t& to) const
  {
    Monomial result;

    BOOST_FOREACH(const typename exponent_map_t::value_type& exponentMapping, exponents)
    {
      if (exponentMapping.first == from)
        result.exponents.insert(std::make_pair(to, exponentMapping.second));
      else
        result.exponents.insert(exponentMapping);
    }

    return result;
  }
};

template<typename V>
std::ostream& operator<<(std::ostream& o, const Monomial<V>& p)
{
  return p.write(o);
}

}

namespace std
{
  template<typename V>
  void swap(cfd::Monomial<V>& a, cfd::Monomial<V>& b)
  {
    a.swap(b);
  }
}

#endif
