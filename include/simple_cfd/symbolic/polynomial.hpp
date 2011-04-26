#ifndef SIMPLE_CFD_SYMBOLIC_POLYNOMIAL_HPP
#define SIMPLE_CFD_SYMBOLIC_POLYNOMIAL_HPP

#include <utility>
#include <algorithm>
#include <vector>
#include <ostream>
#include <cassert>
#include <iterator>
#include <stack>
#include <queue>
#include <functional>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/operators.hpp>
#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/symbol.hpp>

namespace cfd
{

namespace symbolic
{

class Polynomial
{
private:
  class Monomial;
  class MonomialProduct;

  typedef std::pair<Symbol, std::size_t> exponent_t;
  std::vector<exponent_t> exponents;
  std::vector<std::size_t> monomialOffsets;
  std::vector<Rational> coefficients;

  /* This implementation is counter-intuitive. "b" is actually less than
     "a" since they are actually "a^0*b^1" and "a^1", respectively.
  */
  struct ExponentComparator
  {
    bool operator()(const exponent_t& a, const exponent_t& b) const
    {
      if (a.first == b.first)
        return a.second < b.second;
      else
        return a.first > b.first;
    }
  };

  class Monomial : public boost::totally_ordered<Monomial>
  {
  private:
    const exponent_t* beginExponent;
    const exponent_t* endExponent;

  public:
    typedef const exponent_t* iterator;
    typedef const exponent_t* const_iterator;

    Monomial(const exponent_t* _begin, const exponent_t* _end) : 
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
      return end() - begin();
    }

    bool operator<(const Monomial& m) const
    {
      return std::lexicographical_compare(begin(), end(), m.begin(),
        m.end(), ExponentComparator());
    }

    bool operator==(const Monomial& m) const
    {
      return size() == m.size() && std::equal(begin(), end(), m.begin());
    }

    MonomialProduct operator*(const Monomial& m) const;

    void write(std::ostream& out)
    {
      BOOST_FOREACH(const exponent_t& e, *this)
      {
        out << e.first;
        
        if (e.second != 1)
          out << "^{" << e.second << "}";
      }
    }
  };

  class MonomialProductIter : public boost::iterator_facade<MonomialProductIter, const Polynomial::exponent_t, 
                                       std::forward_iterator_tag, const Polynomial::exponent_t>
  {
  private:
    bool smaller;
    Monomial::iterator firstIter, firstEnd;
    Monomial::iterator secondIter, secondEnd;

  public:
    MonomialProductIter(const Monomial& first, const Monomial& second, const bool begin) :
      firstIter(begin ? first.begin() : first.end()), firstEnd(first.end()), 
      secondIter(begin ? second.begin() : second.end()), secondEnd(second.end())
    {
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
      else if (firstIter->first == secondIter->first)
      {
        return reference(firstIter->first, firstIter->second + secondIter->second);
      }
      else
      {
        const bool firstIsSmaller = firstIter->first < secondIter->first;
        return *(firstIsSmaller ? firstIter : secondIter);
      }
    }

    void increment()
    {
      if (firstIter == firstEnd)
      {
        ++secondIter;
      }
      else if (secondIter == secondEnd)
      {
        ++firstIter;
      }
      else if (firstIter->first == secondIter->first)
      {
        ++firstIter;
        ++secondIter;
      }
      else
      {
        const bool firstIsSmaller = firstIter->first < secondIter->first;
        ++(firstIsSmaller ? firstIter : secondIter); 
      }
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
        m.end(), ExponentComparator());
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
  
  class PolyIter : public boost::iterator_facade<PolyIter, 
                            const Monomial, std::random_access_iterator_tag, const Monomial>
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

  const Monomial operator[](const std::size_t index) const
  {
    return Monomial(&exponents[monomialOffsets[index]], 
                    &exponents[monomialOffsets[index+1]]); 
  }

  const Rational coefficient(std::size_t index) const
  {
    return coefficients[index];
  }

  template<typename monomial_t>
  void pushMonomialTemplated(const Rational& c, const monomial_t& m)
  {
    if (numTerms() == 0 || m != (*this)[numTerms()-1])
    {
      coefficients.push_back(c);
      exponents.insert(exponents.end(), m.begin(), m.end());
      monomialOffsets.push_back(monomialOffsets[numTerms()-1] + m.size());
    }
    else
    {
      coefficients[numTerms()-1] += c;
    }
  }

  void pushMonomial(const Rational& c, const Monomial& m)
  {
    pushMonomialTemplated(c, m);
  }

  void pushMonomial(const Rational& c, const MonomialProduct& m)
  {
    pushMonomialTemplated(c, m);
  }
      
public:
  Polynomial()
  {
    monomialOffsets.push_back(0);
  }

  Polynomial(const Rational& r)
  {
    coefficients.push_back(r);
    monomialOffsets.push_back(0);
    monomialOffsets.push_back(0);
  }

  Polynomial(const Symbol& s)
  {
    coefficients.push_back(Rational(1));
    monomialOffsets.push_back(0);
    monomialOffsets.push_back(1);
    exponents.push_back(exponent_t(s, 1));
  }

  Polynomial operator+(const Polynomial& b) const
  {
    Polynomial result;

    const Polynomial& a = *this;
    const_iterator aIter = a.begin();
    const_iterator bIter = b.begin();

    while(aIter != a.end() || bIter != b.end())
    {
      if (aIter == a.end())
      {
        result.pushMonomial(bIter.coefficient(), *bIter);
        ++bIter;
      }
      else if (bIter == b.end())
      {
        result.pushMonomial(aIter.coefficient(), *aIter);
        ++aIter;
      }
      else if (*aIter < *bIter)
      {
        result.pushMonomial(aIter.coefficient(), *aIter);
        ++aIter;
      }
      else
      {
        result.pushMonomial(bIter.coefficient(), *bIter);
        ++bIter;
      }
    }

    return result;
  }

  Polynomial operator*(const Polynomial& b) const
  {
    typedef std::pair<MonomialProduct, std::size_t> state_t;

    Polynomial result;
    std::vector<const_iterator> positions(numTerms(), b.begin());
    std::priority_queue<state_t, 
                        std::vector<state_t>,
                        std::greater<state_t> > queue;

    if (b.begin() != b.end())
    {
      for(std::size_t i=0; i<numTerms(); ++i)
        queue.push(std::make_pair((*this)[i]*(*positions[i]), i));
    }

    while(!queue.empty())
    {
      const std::pair<MonomialProduct, std::size_t> top = queue.top();
      queue.pop();

      const int position = top.second;
      result.pushMonomial(positions[position].coefficient() * coefficients[position], top.first);
      ++positions[position];

      if (positions[position] != b.end())
        queue.push(std::make_pair((*this)[top.second]*(*positions[top.second]), top.second));
    }

    return result;
  }

  std::size_t numTerms() const
  {
    return coefficients.size();
  }

  void write(std::ostream& out) const
  {
    for(const_iterator iter = begin(); iter != end(); ++iter)
    {
      out << iter.coefficient();
      iter->write(out);

      if (iter+1 != end())
        out << " + ";
    }
  }
};

}

}

#endif
