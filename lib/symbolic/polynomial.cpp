#include <queue>
#include <utility>
#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/polynomial.hpp>
#include <simple_cfd/util/hash.hpp>

namespace cfd
{

namespace symbolic
{

Polynomial::Polynomial()
{
  monomialOffsets.push_back(0);
}

Polynomial::Polynomial(const Rational& r)
{
  coefficients.push_back(r);
  monomialOffsets.push_back(0);
  monomialOffsets.push_back(0);
}

Polynomial::Polynomial(const Symbol& s)
{
  coefficients.push_back(Rational(1));
  monomialOffsets.push_back(0);
  monomialOffsets.push_back(1);
  exponents.push_back(exponent_t(s, 1));
}

std::size_t Polynomial::nops() const
{
  return numTerms();
}

bool Polynomial::depends(const std::set<Symbol>& symbols) const
{
  BOOST_FOREACH(const exponent_t& exponent, exponents)
  {
    if (symbols.find(exponent.first) != symbols.end())
      return true;
  }

  return false;
}

std::size_t Polynomial::untypedHash() const
{
  std::size_t result = 0x07b61b4a;

  BOOST_FOREACH(const exponent_t& exponent, exponents)
  {
    util::hash_accum(result, exponent.first);
    util::hash_accum(result, exponent.second);
  }

  BOOST_FOREACH(const std::size_t offset, monomialOffsets)
  {
    util::hash_accum(result, offset);
  }

  BOOST_FOREACH(const Rational& c, coefficients)
  {
    util::hash_accum(result, c);
  }

  return result;
}

Expr Polynomial::derivative(const Symbol& s) const
{
  Polynomial result;
  for(const_iterator iter = begin(); iter != end(); ++iter)
  {
    const long exponent = iter->getExponent(s);

    if (exponent>0)
    {
      const Rational coefficient = iter.coefficient() * exponent;
      const exponent_t decrement(s, -1);
      const Monomial decrementMonomial(&decrement, &decrement + 1);
      result.pushMonomial(coefficient, decrementMonomial * (*iter));
    }
  }

  return result;
}

Expr Polynomial::integrate(const Symbol& s, const unsigned flags) const
{
  Polynomial result;
  for(const_iterator iter = begin(); iter != end(); ++iter)
  {
    const long exponent = iter->getExponent(s);

    const Rational coefficient = iter.coefficient() / (exponent+1);
    const exponent_t increment(s, 1);
    const Monomial incrementMonomial(&increment, &increment + 1);
    result.pushMonomial(coefficient, incrementMonomial * (*iter));
  }

  return result;
}

Float Polynomial::eval(const Expr::subst_map& map) const
{
  Float result = 0.0;
  
  for(const_iterator iter = begin(); iter != end(); ++iter)
  {
    Float term = iter.coefficient().toFloat();
    for(Monomial::const_iterator powerIter = iter->begin(); powerIter != iter->end(); ++powerIter)
    {
      const Float symbolValue = powerIter->first.eval(map);
      const long exponent = powerIter->second;

      term *= pow(symbolValue, exponent);
    }

    result += term;
  }

  return result;
}

void Polynomial::accept(NumericExpressionVisitor<Symbol>& v) const
{
  for(const_iterator iter = begin(); iter != end(); ++iter)
  {
    const Rational coefficient = iter.coefficient();
    const bool hasCoefficient = (coefficient != 1);

    if (hasCoefficient)
      coefficient.accept(v);

    for(Monomial::const_iterator powerIter = iter->begin(); powerIter != iter->end(); ++powerIter)
    {
      v.visitVariable(powerIter->first);

      if (powerIter->second != 1)
        v.visitExponent(powerIter->second);
    }

    v.postProduct(iter->size() + (hasCoefficient ? 1 : 0));
  }

  v.postSummation(numTerms());
}

bool Polynomial::operator==(const Polynomial& b) const
{
  return monomialOffsets == b.monomialOffsets 
         && exponents == b.exponents
         && coefficients == b.coefficients;
}

Polynomial& Polynomial::operator/=(const Rational& r)
{
  this->invalidateHash();

  BOOST_FOREACH(Rational& c, coefficients)
    c /= r;

  return *this;
}

Polynomial& Polynomial::operator*=(const Rational& r)
{
  this->invalidateHash();

  BOOST_FOREACH(Rational& c, coefficients)
    c *= r;

  return *this;
}

Polynomial Polynomial::operator+(const Polynomial& b) const
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

Polynomial::MonomialProduct Polynomial::Monomial::operator*(const Monomial& m) const
{
  return MonomialProduct(*this, m);
}

Polynomial Polynomial::operator*(const Polynomial& b) const
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

Expr Polynomial::toSum() const
{
  Sum result;
  for(const_iterator polyIter = begin(); polyIter != end(); ++polyIter)
  {
    Product term = Product::constant(polyIter.coefficient());

    BOOST_FOREACH(const exponent_t& exponent, *polyIter)
      term *= Product::pow(exponent.first, exponent.second);

    result += term;
  }

  return result.simplify();
}

Expr Polynomial::subs(const Expr::subst_map& map, const unsigned flags) const
{
  return toSum().subs(map, flags);
}

void Polynomial::write(std::ostream& out) const
{
  for(const_iterator iter = begin(); iter != end(); ++iter)
  {
    out << iter.coefficient();
    iter->write(out);

    if (iter+1 != end())
      out << " + ";
  }
}

std::ostream& operator<<(std::ostream& o, const Polynomial& p)
{
  p.write(o);
  return o;
}

}

}
