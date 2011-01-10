#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/number.hpp>
#include <map>
#include <utility>
#include <boost/foreach.hpp>

namespace cfd
{

namespace symbolic
{

Expr Product::null() const
{
  return Expr(new Number(1));
}

void Product::write(std::ostream& o) const
{
  if (begin() == end())
  {
    o << "1.0";
    return;
  }

  const_iterator current(begin());

  o << "(";
  while(current != end())
  {
    const bool expIsOne = current->second == 1;
    o << (expIsOne ? "" : "(") << current->first << (expIsOne ? "" : ")");

    if (!expIsOne) 
      o << "^" << "{" << current->second << "}";

    ++current;

    if (current != end())
      o << "*";
  }
  o << ")";
}

Expr Product::derivative(const Symbol& s) const
{
  Sum summation;

  BOOST_FOREACH(const TermMap::value_type d, std::make_pair(begin(), end()))
  {
    TermMap newTerm(begin(), end());
    newTerm.erase(d.first);
    ++newTerm[Number(d.second)];
    ++newTerm[d.first.derivative(s)];
    newTerm[d.first]+=d.second-1;
    summation = summation + Product(newTerm);
  }
  
  return summation;
}

Expr Product::simplify() const
{
  const Number zero(0);
  const Expr simplified = PairSeq<Product>::simplify();
  const Basic& basic = simplified.internal();

  if (PairSeq<Product>::getType(basic) == PairSeq<Product>::getType(*this))
  {
    const Product& p = static_cast<const Product&>(basic);

    BOOST_FOREACH(const TermMap::value_type termPair, std::make_pair(p.begin(), p.end()))
    {
      if (termPair.first == zero && termPair.second > 0)
        return Number(0);
    }
  }
  
  return simplified;
}

}

}
