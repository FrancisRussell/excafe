#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/number.hpp>
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

}

}
