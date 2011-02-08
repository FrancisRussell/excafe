#ifndef SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP
#define SIMPLE_CFD_SYMBOLIC_PAIR_SEQ_HPP

#include <cstddef>
#include <map>
#include <ostream>
#include <utility>
#include <vector>
#include "abstract_basic.hpp"
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
  typedef std::map<Expr, int> TermMap;

protected:
  typedef T child_type;
  TermMap terms;

protected:
  PairSeq()
  {
  }

  PairSeq(const TermMap& _terms): terms(_terms)
  {
  }

  PairSeq(const Expr& a, const Expr& b) 
  {
    ++terms[a];
    ++terms[b];
  }
  
  virtual Expr null() const = 0;

  void addSimplifiedTerms(const int multiplier, TermMap& newTermMap, const child_type& seq) const
  {
    const Expr nullExpr = null();
    BOOST_FOREACH(const TermMap::value_type term, std::make_pair(seq.begin(), seq.end()))
    {
      // Terms multiplied by 0 or raised to 0 can be removed from the sequence
      if (term.second != 0)
      {
        const Expr simplified = term.first.simplify();

        if (AbstractBasic<T>::getType(simplified.internal()) == AbstractBasic<T>::getType(*this))
        {
          const child_type& child(static_cast<const child_type&>(simplified.internal()));
          addSimplifiedTerms(multiplier*term.second, newTermMap, child);
        }
        else if (simplified != nullExpr)
        {
          // Terms equal that are 0 in a sum, or 1 in a product can be removed
          newTermMap[simplified] += multiplier*term.second;
        }
      }
    }
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

  virtual std::size_t nops() const
  {
    return terms.size();
  }

  bool operator<(const child_type& s) const
  {
    return terms < s.terms;
  }

  bool operator==(const child_type& s) const
  {
    return terms == s.terms;
  }

  Expr simplify() const
  {
    const Expr nullExpr = null();
    TermMap newTermMap;
    addSimplifiedTerms(1, newTermMap, static_cast<const child_type&>(*this));

    if (newTermMap.empty())
    {
      return nullExpr;
    }
    else if (newTermMap.size() == 1 && newTermMap.begin()->second == 1)
    {
      return newTermMap.begin()->first;
    }
    else
    {
      return child_type(newTermMap);
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
    return boost::hash_range(terms.begin(), terms.end());
  }
};

}

}

#endif
