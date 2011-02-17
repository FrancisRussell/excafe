#ifndef SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP
#define SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP

#include <cstddef>
#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>
#include <ostream>
#include <utility>
#include <vector>
#include "abstract_basic.hpp"
#include "rational.hpp"
#include "expr.hpp"
#include <iostream>
#include <boost/foreach.hpp>

namespace cfd
{

namespace symbolic
{

template<typename T>
class PairSeq : public AbstractBasic<T>
{
protected:
  typedef boost::unordered_map<Expr, int> TermMap;

protected:
  typedef T child_type;
  Rational overall;
  TermMap terms;
  bool simplified;

protected:
  PairSeq() : overall(child_type::null()), simplified(false)
  {
  }

  PairSeq(const TermMap& _terms): overall(child_type::null()), terms(_terms), simplified(false)
  {
  }

  PairSeq(const Expr& a, const Expr& b) : overall(child_type::null()), simplified(false)
  {
    ++terms[a];
    ++terms[b];
  }
  
  void addSimplifiedTerms(const int multiplier, TermMap& newTermMap, const child_type& seq) const
  {
    const Expr nullExpr = child_type::null();
    BOOST_FOREACH(const TermMap::value_type term, std::make_pair(seq.begin(), seq.end()))
    {
      const Expr simplified = term.first.simplify();

      if (AbstractBasic<T>::getType(simplified.internal()) == AbstractBasic<T>::getType(*this))
      {
        const child_type& child = static_cast<const child_type&>(simplified.internal());
        newTermMap[child.getOverall()] += multiplier;
        addSimplifiedTerms(multiplier*term.second, newTermMap, child);
      }
      else if (simplified != nullExpr)
      {
        // Terms equal that are 0 in a sum, or 1 in a product can be removed
        newTermMap[simplified] += multiplier*term.second;
      }
    }

    TermMap::iterator newIter(newTermMap.begin());
    while(newIter != newTermMap.end())
    {
      const TermMap::iterator nextIter = boost::next(newIter);

      // Terms multiplied by 0 or raised to 0 can be removed from the sequence
      if (newIter->second == 0)
        newTermMap.erase(newIter);

      newIter = nextIter;
    }
  }

  static void updateOverall(Rational& overall, TermMap& termMap)
  {
    TermMap::iterator iter(termMap.begin());
    while(iter != termMap.end())
    {
      const TermMap::iterator nextIter = boost::next(iter);

      // Terms multiplied by 0 or raised to 0 can be removed from the sequence
      if (is_a<Rational>(iter->first))
      {
        const Rational value = convert_to<Rational>(iter->first);
        const int coefficient = iter->second;
        child_type::combineOverall(overall, child_type::applyCoefficient(value, coefficient));
        termMap.erase(iter);
      }

      iter = nextIter;
    }
  }

  child_type withoutOverall() const
  {
    child_type result(asChild(*this));

    if (overall != child_type::null())
    {
      ++result.terms[result.overall];
      result.overall = child_type::null();
    }
    return result;
  }

public:
  typedef TermMap::value_type value_type;
  typedef TermMap::const_iterator iterator;
  typedef TermMap::const_iterator const_iterator;

  const_iterator begin() const
  {
    return terms.begin();
  }

  const_iterator end() const
  {
    return terms.end();
  }

  Rational getOverall() const
  {
    return overall;
  }

  virtual std::size_t nops() const
  {
    return terms.size();
  }

  bool operator==(const child_type& s) const
  {
    return overall == s.overall
           && terms == s.terms;
  }

  Expr simplify() const
  {
    if (simplified)
      return this->clone();

    TermMap newTermMap;
    Rational newOverall(overall);
    addSimplifiedTerms(1, newTermMap, static_cast<const child_type&>(*this));
    updateOverall(newOverall, newTermMap);

    const Expr nullExpr = child_type::null();
    if (newTermMap.empty())
    {
      return newOverall;
    }
    else if (newTermMap.size() == 1 && newTermMap.begin()->second == 1 && newOverall == child_type::null())
    {
      return newTermMap.begin()->first;
    }
    else
    {
      child_type result(newTermMap);
      result.overall = newOverall;
      result.simplified = true;
      return result;
    }
  }

  Expr subs(const Expr::subst_map& map) const
  {
    TermMap newTermMap;
    BOOST_FOREACH(const TermMap::value_type term, std::make_pair(begin(), end()))
    {
      newTermMap[term.first.subs(map)] += term.second;
    }
    return child_type(newTermMap);
  }

  bool has(const Expr& e) const
  {
    if (e == *this)
      return true;

    BOOST_FOREACH(const TermMap::value_type term, std::make_pair(begin(), end()))
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
    boost::hash_range(result, terms.begin(), terms.end());
    return result;
  }
};

}

}

#endif
