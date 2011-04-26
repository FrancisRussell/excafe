#include <queue>
#include <utility>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/polynomial.hpp>

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

Expr Polynomial::toExpr() const
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
