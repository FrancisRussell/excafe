#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/number.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <map>
#include <utility>
#include <boost/foreach.hpp>

namespace cfd
{

namespace symbolic
{

Expr Sum::null() const
{
  return Expr(new Number(0));
}

void Sum::write(std::ostream& o) const
{
  if (begin() == end())
  {
    o << "0.0";
    return;
  }

  const_iterator current(begin());

  o << "(";
  while(current != end())
  {
    o << (current->second < 0 ? "-" : "");

    if (std::abs(current->second) != 1)
      o << std::abs(current->second) << "*";

    o << current->first;
    ++current;

    if (current != end())
      o << " + ";
  }
  o << ")";
}

Sum Sum::operator+(const Expr& e) const
{
  TermMap map(begin(), end());
  ++map[e];
  return Sum(map);
}

Expr Sum::derivative(const Symbol& s) const
{
  TermMap newTerms;
  BOOST_FOREACH(const TermMap::value_type& e, std::make_pair(begin(), end()))
  {
    newTerms.insert(std::make_pair(e.first.derivative(s), e.second));
  }
  return Sum(newTerms);
}

Expr Sum::integrate(const Symbol& s) const
{
  TermMap newTerms;
  BOOST_FOREACH(const TermMap::value_type& e, std::make_pair(begin(), end()))
  {
    newTerms.insert(std::make_pair(e.first.integrate(s), e.second));
  }
  return Sum(newTerms);
}

Sum Sum::expandedProduct(const Sum& s) const
{
  TermMap newTerms;

  BOOST_FOREACH(const TermMap::value_type& a, std::make_pair(begin(), end()))
  {
    BOOST_FOREACH(const TermMap::value_type& b, std::make_pair(s.begin(), s.end()))
    {
      newTerms[Product(a.first, b.first)] += a.second*b.second;
    }
  }

  return Sum(newTerms);
}

}

}
