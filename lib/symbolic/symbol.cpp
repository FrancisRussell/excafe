#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <boost/functional/hash.hpp>
#include <iostream>

namespace cfd
{

namespace symbolic
{

int Symbol::nextSerial = 0;

Symbol::Symbol(const std::string& _name) : name(_name), serial(nextSerial++)
{
}

std::size_t Symbol::nops() const
{
  return 0;
}

void Symbol::write(std::ostream& o) const
{
  o << name;
}

Expr Symbol::derivative(const Symbol& s) const
{
  if (serial == s.serial)
    return Expr(new Rational(1));
  else
    return Expr(new Rational(0));
}

bool Symbol::operator==(const Symbol& s) const
{
  return serial == s.serial;
}

bool Symbol::operator<(const Symbol& s) const
{
  return serial < s.serial;
}

bool Symbol::has(const Expr& e) const
{
  return e == *this;
}

std::size_t Symbol::untypedHash() const
{
  return boost::hash<int>()(serial);
}

Expr Symbol::subs(const Expr::subst_map& map) const
{
  const Expr::subst_map::const_iterator iter = map.find(*this);

  if (iter != map.end())
  {
    return iter->second;
  }
  else
  {
    return clone();
  }
}

Expr Symbol::eval() const
{
  return clone();
}

Expr Symbol::integrate(const Symbol& s) const
{
  if (serial != s.serial)
    return Product(*this, s);
  else
    return Product(Rational(1,2), Product::pow(s, 2));
}

void Symbol::accept(NumericExpressionVisitor<Symbol>& v) const
{
  v.visitVariable(*this);
}

}

}