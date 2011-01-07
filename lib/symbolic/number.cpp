#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/number.hpp>
#include <boost/functional/hash.hpp>

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

Expr Number::clone() const
{
  return Expr(new Number(*this));
}

Expr Number::derivative(const Symbol& s) const
{
  return Expr(new Number(0));
}

bool Number::isNumber() const
{
  return true;
}

std::size_t Number::hashValue() const
{
  return boost::hash<double>()(value);
}


}

}
