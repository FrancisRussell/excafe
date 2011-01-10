#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/number.hpp>
#include <boost/functional/hash.hpp>

#include <iostream>
namespace cfd
{

namespace symbolic
{

Number::Number(const double _value) : value(_value)
{
}

std::size_t Number::nops() const
{
  return 0;
}

void Number::write(std::ostream& o) const
{
  o << value;
}

Expr Number::derivative(const Symbol& s) const
{
  return Expr(new Number(0));
}

bool Number::isNumber() const
{
  return true;
}

bool Number::operator<(const Number& n) const
{
  return value < n.value;
}

bool Number::operator==(const Number& n) const
{
  return value == n.value;
}

std::size_t Number::untypedHash() const
{
  return boost::hash<double>()(value);
}

Expr Number::subs(const Expr::subst_map& map) const
{
  return clone();
}

}

}
