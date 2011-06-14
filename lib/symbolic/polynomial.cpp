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
  if (r != 0)
  {
    coefficients.push_back(r);
    monomialOffsets.push_back(0);
  }
  monomialOffsets.push_back(0);
}

Polynomial::Polynomial(const Symbol& s)
{
  coefficients.push_back(Rational(1));
  monomialOffsets.push_back(0);
  monomialOffsets.push_back(1);
  exponents.push_back(exponent_t(s, 1));
}

Expr Polynomial::simplify() const
{
  if (monomialOffsets.size() == 2)
  {
    if (monomialOffsets[1] == 0)
      return coefficients[0];

    if (monomialOffsets[1] == 1
        && coefficients[0] == 1
        && exponents[0].exp() == 1)
      return exponents[0].sym();
  }
  
  return toSum();
}

std::size_t Polynomial::nops() const
{
  return numTerms();
}

bool Polynomial::depends(const std::set<Symbol>& symbols) const
{
  BOOST_FOREACH(const exponent_t& exponent, exponents)
  {
    if (symbols.find(exponent.sym()) != symbols.end())
      return true;
  }

  return false;
}

std::size_t Polynomial::untypedHash() const
{
  std::size_t result = 0x07b61b4a;

  BOOST_FOREACH(const exponent_t& exponent, exponents)
  {
    util::hash_accum(result, exponent.sym());
    util::hash_accum(result, exponent.exp());
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

Expr Polynomial::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  Polynomial result;
  for(const_iterator iter = begin(); iter != end(); ++iter)
  {
    const long exponent = iter->getExponent(s);

    if (exponent>0)
    {
      const Rational coefficient = iter.coefficient() * exponent;
      const exponent_t decrement(s, -1);
      const detail::Monomial decrementMonomial(&decrement, &decrement + 1);
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
    const detail::Monomial incrementMonomial(&increment, &increment + 1);
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
    for(detail::Monomial::const_iterator powerIter = iter->begin(); powerIter != iter->end(); ++powerIter)
    {
      const Float symbolValue = powerIter->sym().eval(map);
      const long exponent = powerIter->exp();

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

    for(detail::Monomial::const_iterator powerIter = iter->begin(); powerIter != iter->end(); ++powerIter)
    {
      v.visitVariable(powerIter->sym());

      if (powerIter->exp() != 1)
        v.visitExponent(powerIter->exp());
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
  if (r == 1)
    return *this;

  this->invalidateHash();

  BOOST_FOREACH(Rational& c, coefficients)
    c /= r;

  return *this;
}

Rational Polynomial::findMultiplier() const
{
  Rational result = 0;

  BOOST_FOREACH(const Rational& c, coefficients)
  {
    result = Rational::gcd(result, c);

    if (result == 1)
      break;
  }

  return result;
}

Polynomial& Polynomial::operator*=(const Rational& r)
{
  if (r == 1)
    return *this;

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

detail::MonomialProduct detail::Monomial::operator*(const detail::Monomial& m) const
{
  return MonomialProduct(*this, m);
}

Polynomial Polynomial::operator*(const Polynomial& b) const
{
  typedef std::pair<const_iterator, const_iterator> state_t;

  Polynomial result;
  std::priority_queue<state_t, 
                      std::vector<state_t>,
                      MultiplyComparator> queue;

  if (b.begin() != b.end())
  {
    for(const_iterator iter=begin(); iter != end(); ++iter)
      queue.push(state_t(iter, b.begin()));
  }

  while(!queue.empty())
  {
    state_t top = queue.top();
    queue.pop();

    result.pushMonomial(top.first.coefficient() * top.second.coefficient(), *top.first * *top.second);
    ++top.second;

    if (top.second != b.end())
      queue.push(top);
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
      term *= Product::pow(exponent.sym(), exponent.exp());

    result += term;
  }

  return result.simplify();
}

Expr Polynomial::subs(const Expr::subst_map& map, const unsigned flags) const
{
  bool hasNonRationals = false;

  BOOST_FOREACH(const Expr::subst_map::value_type& sub, map)
  {
    if (!is_exactly_a<Rational>(sub.second))
    {
      hasNonRationals = true;
      break;
    }
  }

  if (hasNonRationals)
  {
    return toSum().subs(map, flags);
  }
  else
  {
    Polynomial result;
    for(const_iterator polyIter = begin(); polyIter != end(); ++polyIter)
    {
      Rational newCoeff = polyIter.coefficient();
      std::vector<exponent_t> exponents;

      for(detail::Monomial::const_iterator powerIter = polyIter->begin(); powerIter != polyIter->end(); ++powerIter)
      {
        const Expr::subst_map::const_iterator subIter = map.find(powerIter->sym());
        if (subIter != map.end())
        {
          const Rational& r = convert_to<Rational>(subIter->second);
          newCoeff *= pow(r, powerIter->exp());
        }
        else
        {
          exponents.push_back(*powerIter);
        }
      }

      const detail::Monomial newMonomial(&(*exponents.begin()), &(*exponents.end()));
      result.pushMonomial(newCoeff, newMonomial);
    }
    return result;
  }
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

Expr Polynomial::extractPolynomials(ExtractedExpressions& extracted) const
{
  return clone();
}

std::ostream& operator<<(std::ostream& o, const Polynomial& p)
{
  p.write(o);
  return o;
}

}

}
