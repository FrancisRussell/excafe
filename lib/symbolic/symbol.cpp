#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/util/hash.hpp>
#include <simple_cfd/exception.hpp>
#include <set>

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

bool Symbol::depends(const std::set<Symbol>& symbols) const
{
  return symbols.find(*this) != symbols.end();
}

std::size_t Symbol::untypedHash() const
{
  std::size_t result = 0x2d3a117b;
  cfd::util::hash_accum(result, serial);
  return result;
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

Float Symbol::eval(const Expr::subst_map& map) const
{
  const Expr::subst_map::const_iterator iter = map.find(*this);

  if (iter != map.end() && is_exactly_a<Float>(iter->second))
  {
    return convert_to<Float>(iter->second);
  }
  else
  {
    CFD_EXCEPTION("Missing value in evaluation map.");
  }
}

Expr Symbol::integrate_internal(const Symbol& s) const
{
  if (serial != s.serial)
    return Product::mul(*this, s);
  else
    return Product::mul(Rational(1,2), Product::pow(s, 2));
}

void Symbol::accept(NumericExpressionVisitor<Symbol>& v) const
{
  v.visitVariable(*this);
}

}

}
