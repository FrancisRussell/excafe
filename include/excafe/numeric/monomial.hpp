#ifndef EXCAFE_NUMERIC_MONOMIAL_HPP
#define EXCAFE_NUMERIC_MONOMIAL_HPP

#include <map>
#include <cstddef>
#include <cmath>
#include <set>
#include <iosfwd>
#include <utility>
#include <ostream>
#include <boost/foreach.hpp>
#include <boost/operators.hpp>
#include <excafe/util/lazy_copy.hpp>
#include "value_map.hpp"
#include "cast.hpp"

namespace excafe
{

template<typename V, typename T>
class Monomial : boost::multipliable<Monomial<V,T>,
                 boost::totally_ordered<Monomial<V,T>
                 > >
{
public:
  typedef V variable_t;
  typedef T numeric_type;

private:
  typedef std::map<variable_t, std::size_t> exponent_map_t;
  util::LazyCopy<exponent_map_t> exponents;

public:
  typedef typename exponent_map_t::value_type value_type;
  typedef typename exponent_map_t::const_iterator iterator;
  typedef typename exponent_map_t::const_iterator const_iterator;

  Monomial()
  {
  }

  Monomial(const Monomial& m) : exponents(m.exponents)
  {
  }

  Monomial(const variable_t& variable, const std::size_t exponent)
  {
    if (exponent != 0)
      exponents->insert(std::make_pair(variable, exponent));
  }

  const_iterator begin() const 
  {
    return exponents->begin();
  }

  const_iterator end() const 
  {
    return exponents->end();
  }

  std::size_t size() const
  {
    return exponents->size();
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

  std::pair<numeric_type, Monomial> derivative(const variable_t& variable) const
  {
    const numeric_type coefficient = excafe::numeric_cast<double>(getExponent(variable));
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

  std::pair<numeric_type, Monomial> substituteValues(const detail::ValueMap<variable_t, numeric_type>& valueMap) const
  {
    Monomial result;
    numeric_type coefficient = 1.0;

    typedef typename detail::ValueMap<variable_t, numeric_type>::var_subst_map    var_subst_map;
    typedef typename detail::ValueMap<variable_t, numeric_type>::scalar_subst_map scalar_subst_map;

    const var_subst_map&    varSubstMap    = valueMap.getVariableSubstitutions();
    const scalar_subst_map& scalarSubstMap = valueMap.getScalarSubstitutions();

    BOOST_FOREACH(const typename exponent_map_t::value_type& exponentMapping, *exponents)
    {
      const typename scalar_subst_map::const_iterator scalarSubstIter = scalarSubstMap.find(exponentMapping.first);
      const typename var_subst_map::const_iterator varSubstIter = varSubstMap.find(exponentMapping.first);

      if (scalarSubstIter != scalarSubstMap.end())
      {
        coefficient *= pow(scalarSubstIter->second, exponentMapping.second);
      }
      else if (varSubstIter != varSubstMap.end())
      {
        (*result.exponents)[varSubstIter->second] += exponentMapping.second;
      }
      else
      {
        (*result.exponents)[exponentMapping.first] += exponentMapping.second;
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

    BOOST_FOREACH(const typename exponent_map_t::value_type& exponentMapping, *exponents)
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
  void swap(excafe::Monomial<V,T>& a, excafe::Monomial<V,T>& b)
  {
    a.swap(b);
  }
}

#endif
