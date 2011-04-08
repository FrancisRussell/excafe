#ifndef SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP
#define SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>
#include <set>
#include <boost/foreach.hpp>
#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>
#include "abstract_basic.hpp"
#include "rational.hpp"
#include "expr.hpp"
#include <simple_cfd/util/lazy_copy.hpp>
#include <simple_cfd/util/hash.hpp>
#include <simple_cfd/util/hash_unordered_map.hpp>

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

  enum RewriteState
  {
    NON_NORMALISED,
    NORMALISED,
    NORMALISED_AND_EXTRACTED
  };

  Rational overall;
  LazyTermMap terms;
  RewriteState rewriteState;

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
  static Expr constructSimplifiedExpr(const Rational& overall, const LazyTermMap& map, 
                                      const RewriteState rewriteState)
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
      result.rewriteState = rewriteState;
      return result;
    }
  }

  PairSeq() : overall(child_type::null()), rewriteState(NON_NORMALISED)
  {
  }

  PairSeq(const Rational& _overall, const LazyTermMap& _terms) : 
    overall(_overall), terms(_terms), rewriteState(NON_NORMALISED)
  {
    if (hasZeros(terms.cref()))
      removeZeros(getTerms());
  }

  RewriteState getRewriteState() const
  {
    return rewriteState;
  }

  void mergeSubTerms(Rational& overall, TermMap& newTermMap, const coeff_type& multiplier, const child_type& seq) const
  {
    child_type::combineOverall(overall, child_type::applyCoefficient(seq.overall, multiplier));

    BOOST_FOREACH(const typename TermMap::value_type& term, seq)
    {
      const Expr localTerm = term.first;
      const coeff_type localMultiplier = multiplier*term.second;

      if (is_exactly_a<child_type>(localTerm))
      {
        const child_type& child = convert_to<child_type>(localTerm.internal());
        mergeSubTerms(overall, newTermMap, localMultiplier, child);
      }
      else if (is_exactly_a<Rational>(localTerm))
      {
        const Rational value = convert_to<Rational>(localTerm);
        child_type::combineOverall(overall, child_type::applyCoefficient(value, localMultiplier));
      }
      else
      {
        newTermMap[localTerm] += localMultiplier;
      }
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

  child_type getNormalised() const
  {
    // The null co-efficient for products and sums. This is 1 for both
    // products and sums, unlike the null value, which is 0 for sums.
    const coeff_type defaultCoefficient(1);
    
    // Ordering of these transformations is important. We extract
    // multipliers first since this causes (possibly nested) single-term
    // sums and products to be simplified to their singleton term. If we
    // don't do this first, mergeSubTerms may miss incorporating
    // child terms.

    const child_type simplifiedChildren = asChild(*this).extractMultipliers();
    LazyTermMap newTermMap;
    Rational newOverall = child_type::null();
    mergeSubTerms(newOverall, *newTermMap, defaultCoefficient, simplifiedChildren);
    removeZeros(*newTermMap);

    return child_type(newOverall, newTermMap);
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
    if (getRewriteState() == NORMALISED || getRewriteState() == NORMALISED_AND_EXTRACTED)
      return this->clone();

    const child_type normalised = getNormalised();
    return constructSimplifiedExpr(normalised.overall, normalised.terms, NORMALISED);
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

  bool depends(const std::set<Symbol>& symbols) const
  {
    BOOST_FOREACH(const typename TermMap::value_type& term, std::make_pair(begin(), end()))
    {
      if (term.first.depends(symbols))
        return true;
    }

    return false;
  }

  std::size_t untypedHash() const
  {
    std::size_t result = 0x7730fe1a;
    cfd::util::hash_accum(result, overall);
    cfd::util::hash_accum(result, cfd::util::hash_unordered_map(getTerms()));
    return result;
  }
};

}

}

#endif
