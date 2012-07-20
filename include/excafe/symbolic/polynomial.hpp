#ifndef EXCAFE_SYMBOLIC_POLYNOMIAL_HPP
#define EXCAFE_SYMBOLIC_POLYNOMIAL_HPP

#include <utility>
#include <algorithm>
#include <vector>
#include <set>
#include <ostream>
#include <cassert>
#include <iterator>
#include <functional>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include "symbol.hpp"
#include "expr.hpp"
#include "abstract_basic.hpp"

namespace excafe
{

namespace symbolic
{


namespace detail
{

template<bool is_value> class Exponent;
class Monomial;
class MonomialProduct;
class MonomialProductIter;

template<bool is_value>
class Exponent
{
private:
  typedef typename boost::mpl::if_<boost::mpl::bool_<is_value>,
                                   Symbol,
                                   const Symbol&>::type symbol_t;

  symbol_t symbol;
  long exponent;

public:
  Exponent(const Symbol& _symbol, const long _exponent) :
    symbol(_symbol), exponent(_exponent)
  {
  }

  Exponent(const Exponent<true>& value) :
    symbol(value.sym()), exponent(value.exp())
  {
  }

  Exponent(const Exponent<false>& ref) :
    symbol(ref.sym()), exponent(ref.exp())
  {
  }

  bool operator==(const Exponent<true>& value) const
  {
    return exponent == value.exp() &&
           symbol == value.sym();
  }

  bool operator==(const Exponent<false>& ref) const
  {
    return exponent == ref.exp() &&
           symbol == ref.sym();
  }

  const Symbol& sym() const
  {
    return symbol;
  }

  const long exp() const
  {
    return exponent;
  }
};

/* This implementation is counter-intuitive. "b" is actually less than
   "a" since they are actually "a^0*b^1" and "a^1", respectively.
*/
struct ExponentComparator
{
  template<bool a_is_value, bool b_is_value>
  bool operator()(const Exponent<a_is_value>& a, const Exponent<b_is_value>& b) const
  {
    if (a.sym() == b.sym())
      return a.exp() < b.exp();
    else
      return a.sym() > b.sym();
  }
};

class Monomial : public boost::totally_ordered<Monomial>
{
private:
  const Exponent<true>* beginExponent;
  const Exponent<true>* endExponent;

public:
  typedef const Exponent<true>* iterator;
  typedef const Exponent<true>* const_iterator;

  Monomial(const Exponent<true>* _begin, const Exponent<true>* _end) : 
    beginExponent(_begin), endExponent(_end)
  {
  }

  const_iterator begin() const
  {
    return beginExponent;
  }

  const_iterator end() const
  {
    return endExponent;
  }

  std::size_t size() const
  {
    return endExponent - beginExponent;
  }

  long getExponent(const Symbol& s) const
  {
    const_iterator iter = begin();

    while(iter != end() && iter->sym() != s)
      ++iter;

    if (iter != end())
      return iter->exp();
    else
      return 0;
  }

  bool operator<(const Monomial& m) const
  {
    return std::lexicographical_compare(begin(), end(), m.begin(),
      m.end(), detail::ExponentComparator());
  }

  bool operator==(const Monomial& m) const
  {
    return size() == m.size() && std::equal(begin(), end(), m.begin());
  }

  MonomialProduct operator*(const Monomial& m) const;

  void write(std::ostream& out) const
  {
    BOOST_FOREACH(const Exponent<true>& e, *this)
    {
      out << e.sym();
      
      if (e.exp() != 1)
        out << "^{" << e.exp() << "}";
    }
  }
};

class MonomialProductIter : public boost::iterator_facade<MonomialProductIter, Exponent<true>, 
                                     std::forward_iterator_tag, Exponent<false> >
{
private:
  Monomial::iterator firstIter, firstEnd;
  Monomial::iterator secondIter, secondEnd;

  void incrementInternal()
  {
    if (firstIter == firstEnd)
    {
      ++secondIter;
    }
    else if (secondIter == secondEnd)
    {
      ++firstIter;
    }
    else if (firstIter->sym() == secondIter->sym())
    {
      ++firstIter;
      ++secondIter;
    }
    else
    {
      const bool firstIsSmaller = firstIter->sym() < secondIter->sym();
      ++(firstIsSmaller ? firstIter : secondIter); 
    }
  }

  void skipZeros()
  {
    while(firstIter != firstEnd 
          && secondIter != secondEnd
          && firstIter->sym() == secondIter->sym() 
          && firstIter->exp() + secondIter->exp() == 0)
    {
      ++firstIter;
      ++secondIter;
    }
  }

public:
  MonomialProductIter(const Monomial& first, const Monomial& second, const bool begin) :
    firstIter(begin ? first.begin() : first.end()), firstEnd(first.end()), 
    secondIter(begin ? second.begin() : second.end()), secondEnd(second.end())
  {
    skipZeros();
  }

  reference dereference() const
  {
    if (firstIter == firstEnd)
    {
      return *secondIter;
    }
    else if (secondIter == secondEnd)
    {
      return *firstIter;
    }
    else if (firstIter->sym() == secondIter->sym())
    {
      return reference(firstIter->sym(), firstIter->exp() + secondIter->exp());
    }
    else
    {
      const bool firstIsSmaller = firstIter->sym() < secondIter->sym();
      return *(firstIsSmaller ? firstIter : secondIter);
    }
  }

  void increment()
  {
    incrementInternal();
    skipZeros();
  }

  bool equal(const MonomialProductIter& i) const
  {
    return firstIter == i.firstIter && secondIter == i.secondIter;
  }
};

class MonomialProduct : public boost::totally_ordered<MonomialProduct>,
                        public boost::totally_ordered<MonomialProduct, Monomial>
{
private:
  Monomial first;
  Monomial second;

  template<typename monomial_t>
  bool equal(const monomial_t& m) const
  {
    return size() == m.size()
           && std::equal(begin(), end(), m.begin());
  }

  template<typename monomial_t>
  bool less(const monomial_t& m) const
  {
    return std::lexicographical_compare(begin(), end(), m.begin(),
      m.end(), detail::ExponentComparator());
  }

public:
  typedef MonomialProductIter iterator;
  typedef MonomialProductIter const_iterator;

  MonomialProduct(const Monomial& _first, const Monomial& _second) : first(_first), second(_second)
  {
  }

  iterator begin() const
  {
    return iterator(first, second, true);
  }

  iterator end() const
  {
    return iterator(first, second, false);
  }

  std::size_t size() const
  {
    return std::distance(begin(), end());
  }

  bool operator==(const MonomialProduct& p) const
  {
    return this->equal(p);
  }

  bool operator<(const MonomialProduct& p) const
  {
    return this->less(p);
  }

  bool operator==(const Monomial& m) const
  {
    return this->equal(m);
  }

  bool operator<(const Monomial& m) const
  {
    return this->less(m);
  }
};

}

class Polynomial : public AbstractBasic<Polynomial>,
                   boost::equality_comparable<Polynomial>
{
private:
  typedef detail::Exponent<true> exponent_t;
  std::vector<exponent_t> exponents;
  std::vector<std::size_t> monomialOffsets;
  std::vector<Rational> coefficients;
 
  class PolyIter : public boost::iterator_facade<PolyIter, 
                          detail::Monomial, std::random_access_iterator_tag, const detail::Monomial>
  {
  private:
    const Polynomial* polynomial;
    std::size_t position;

  public:
    PolyIter(const Polynomial& _polynomial, const std::size_t _position) :
      polynomial(&_polynomial), position(_position)
    {
    }

    reference dereference() const
    {
      return (*polynomial)[position];
    }

    Rational coefficient() const
    {
      return polynomial->coefficient(position);
    }

    void increment()
    {
      ++position;
    }

    void decrement()
    {
      --position;
    }

    void advance(const std::size_t n)
    {
      position += n;
    }

    difference_type distance_to(const PolyIter& i) const
    {
      return static_cast<difference_type>(i.position) - position;
    }

    bool equal(const PolyIter& i) const
    {
      return polynomial == i.polynomial && position == i.position;
    }
  };

  typedef PolyIter iterator;
  typedef PolyIter const_iterator;

  const_iterator begin() const
  {
    return PolyIter(*this, 0);
  }

  const_iterator end() const
  {
    return PolyIter(*this, numTerms());
  }

  const detail::Monomial operator[](const std::size_t index) const
  {
    return detail::Monomial(&exponents[monomialOffsets[index]], 
                            &exponents[monomialOffsets[index+1]]); 
  }

  const Rational coefficient(std::size_t index) const
  {
    return coefficients[index];
  }

  template<typename monomial_t>
  void pushMonomialTemplated(const Rational& c, const monomial_t& m)
  {
    if (c == 0)
      return;

    if (numTerms() == 0 || m != (*this)[numTerms()-1])
    {
      coefficients.push_back(c);
      exponents.insert(exponents.end(), m.begin(), m.end());
      monomialOffsets.push_back(monomialOffsets[numTerms()-1] + m.size());
    }
    else
    {
      coefficients[numTerms()-1] += c;

      // Detect cancellation of coefficients.
      if (coefficients[numTerms()-1] == 0)
      {
        exponents.erase(exponents.begin() + monomialOffsets[numTerms()-1], exponents.end());
        coefficients.resize(coefficients.size()-1);
        monomialOffsets.resize(monomialOffsets.size()-1);
      }
    }
  }

  void pushMonomial(const Rational& c, const detail::Monomial& m)
  {
    pushMonomialTemplated(c, m);
  }

  void pushMonomial(const Rational& c, const detail::MonomialProduct& m)
  {
    pushMonomialTemplated(c, m);
  }

  struct MultiplyComparator
  {
    bool operator()(const std::pair<const_iterator, const_iterator>& a,
                    const std::pair<const_iterator, const_iterator>& b) const
    {
      if (a.first >= b.first && a.second >= b.second && a != b)
        return true;
  
      return (*a.first * *a.second) > (*b.first * *b.second);
    }
  };
      
public:
  Polynomial();
  Polynomial(const Rational& r);
  Polynomial(const Symbol& s);

  std::size_t numTerms() const
  {
    return coefficients.size();
  }

  Rational findMultiplier() const;
  Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const;
  Expr integrate(const Symbol& s, unsigned flags) const;
  Expr simplify() const;
  Float eval(const Expr::subst_map& map) const;
  Expr subs(const Expr::subst_map& map, unsigned flags) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  std::size_t nops() const;
  std::size_t untypedHash() const;
  bool depends(const std::set<Symbol>& symbols) const;
  bool operator==(const Polynomial& b) const;
  Polynomial operator+(const Polynomial& b) const;
  Polynomial operator*(const Polynomial& b) const;
  Polynomial& operator/=(const Rational& r);
  Polynomial& operator*=(const Rational& r);
  void write(std::ostream& out) const;
  Expr toSum() const;
  bool isPolynomial() const;
  Expr extractPolynomials(ExtractedExpressions& extracted) const;
};

std::ostream& operator<<(std::ostream& o, const Polynomial& p);

}

}

#endif
