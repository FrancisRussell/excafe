#ifndef SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP
#define SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>
#include "abstract_basic.hpp"
#include "rational.hpp"
#include "expr.hpp"
#include <simple_cfd/util/lazy_copy.hpp>

namespace cfd
{

namespace symbolic
{

template<typename T, typename C>
class PairSeq : public AbstractBasic<T>
{
protected:
  typedef T child_type;
  typedef C coeff_type;
  typedef boost::unordered_map<Expr, coeff_type> TermMap;
  typedef util::LazyCopy<TermMap> LazyTermMap;

  Rational overall;
  LazyTermMap terms;
  bool simplified;

private:
  static bool hasZeros(const TermMap& map)
  {
    BOOST_FOREACH(const typename TermMap::value_type& term, map)
    {
      if (term.second == 0)
        return true;
    }

    return false;
  }

  static void removeZeros(TermMap& map)
  {
    typename TermMap::iterator iter(map.begin());
    while(iter != map.end())
    {
      const typename TermMap::iterator nextIter = boost::next(iter);

      // Terms multiplied by 0 or raised to 0 can be removed from the sequence
      if (iter->second == 0)
        map.erase(iter);

      iter = nextIter;
    }
  }

protected:
  static Expr constructSimplifiedExpr(const Rational& overall, const LazyTermMap& map)
  {
    if (map->empty())
    {
      return overall;
    }
    else if (map->size() == 1 && map->begin()->second == 1 && overall == child_type::null())
    {
      return map->begin()->first;
    }
    else
    {
      child_type result(overall, map);
      result.simplified = true;
      return result;
    }
  }

  PairSeq() : overall(child_type::null()), simplified(false)
  {
  }

  PairSeq(const Rational& _overall, const TermMap& _terms): overall(_overall), terms(_terms), simplified(false)
  {
    if (hasZeros(terms.cref()))
      removeZeros(getTerms());
  }

  PairSeq(const Rational& _overall, const LazyTermMap& _terms): overall(_overall), terms(_terms), simplified(false)
  {
    if (hasZeros(terms.cref()))
      removeZeros(getTerms());
  }

  void addSimplifiedTerms(Rational& overall, TermMap& newTermMap, const coeff_type& multiplier, const child_type& seq) const
  {
    const Expr nullExpr = child_type::null();

    child_type::combineOverall(overall, child_type::applyCoefficient(seq.overall, multiplier));
    BOOST_FOREACH(const typename TermMap::value_type& term, seq)
    {
      const Expr simplified = term.first.simplify();

      if (is_a<child_type>(simplified.internal()))
      {
        const child_type& child = convert_to<child_type>(simplified.internal());
        addSimplifiedTerms(overall, newTermMap, multiplier*term.second, child);
      }
      else if (simplified != nullExpr)
      {
        // Terms equal that are 0 in a sum, or 1 in a product can be removed
        newTermMap[simplified] += multiplier*term.second;
      }
    }
  }

  static void updateOverall(Rational& overall, TermMap& termMap)
  {
    typename TermMap::iterator iter(termMap.begin());
    while(iter != termMap.end())
    {
      const typename TermMap::iterator nextIter = boost::next(iter);

      // Terms multiplied by 0 or raised to 0 can be removed from the sequence
      if (is_a<Rational>(iter->first))
      {
        const Rational value = convert_to<Rational>(iter->first);
        const coeff_type coefficient = iter->second;
        child_type::combineOverall(overall, child_type::applyCoefficient(value, coefficient));
        termMap.erase(iter);
      }

      iter = nextIter;
    }
  }

  void combine(const child_type& other)
  {
    this->invalidateHash();

    if (this == &other)
    {
      child_type::combineOverall(overall, overall);
      BOOST_FOREACH(typename TermMap::value_type& term, getTerms())
      {
        term.second *= 2;
      }
    }
    else
    {
      child_type::combineOverall(overall, other.overall);
      BOOST_FOREACH(const typename TermMap::value_type& term, other.getTerms())
      {
        getTerms()[term.first] += term.second;
      }
    }
  }

  TermMap& getTerms()
  {
    return *terms;
  }

  const TermMap& getTerms() const
  {
    return *terms;
  }

public:
  typedef typename TermMap::value_type value_type;
  typedef typename TermMap::const_iterator iterator;
  typedef typename TermMap::const_iterator const_iterator;

  const_iterator begin() const
  {
    return getTerms().begin();
  }

  const_iterator end() const
  {
    return getTerms().end();
  }

  Rational getOverall() const
  {
    return overall;
  }

  virtual std::size_t nops() const
  {
    return getTerms().size();
  }

  bool operator==(const child_type& s) const
  {
    return overall == s.overall
           && getTerms() == s.getTerms();
  }

  child_type withoutOverall() const
  {
    child_type result(asChild(*this));

    if (overall != child_type::null())
    {
      ++result.getTerms()[result.overall];
      result.overall = child_type::null();
    }
    return result;
  }

  Expr simplify() const
  {
    if (simplified)
      return this->clone();

    // The null co-efficient for products and sums. This is 1 for both
    // products and sums, unlike the null value, which is 0 for sums.
    const coeff_type defaultCoefficient(1);

    LazyTermMap newTermMap;
    Rational newOverall = child_type::null();
    addSimplifiedTerms(newOverall, *newTermMap, defaultCoefficient, asChild(*this));
    updateOverall(newOverall, *newTermMap);
    child_type::extractMultipliers(newOverall, *newTermMap);
    removeZeros(*newTermMap);

    return constructSimplifiedExpr(newOverall, newTermMap);
  }

  Expr subs(const Expr::subst_map& map) const
  {
    LazyTermMap newTermMap;
    BOOST_FOREACH(const typename TermMap::value_type& term, std::make_pair(begin(), end()))
    {
      (*newTermMap)[term.first.subs(map)] += term.second;
    }
    return child_type(overall, newTermMap);
  }

  bool has(const Expr& e) const
  {
    if (e == *this)
      return true;

    BOOST_FOREACH(const typename TermMap::value_type& term, std::make_pair(begin(), end()))
    {
      if (term.first.has(e))
        return true;
    }

    return false;
  }

  std::size_t untypedHash() const
  {
    std::size_t result = 0;
    boost::hash_combine(result, overall);
    boost::hash_range(result, getTerms().begin(), getTerms().end());
    return result;
  }
};

}

}

#endif
