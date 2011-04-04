#include <simple_cfd/symbolic/group.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <ostream>
#include <cstddef>

namespace cfd
{

namespace symbolic
{

Group::Group(const Expr& e) : expr(e)
{
}

std::size_t Group::nops() const
{
  return 1;
}

void Group::write(std::ostream& o) const
{
  o << "group(" << expr << ")";
}

Expr Group::derivative(const Symbol& s) const
{
  return Group(expr.derivative(s)).clone();
}

bool Group::has(const Expr& e) const
{
  return expr.has(e);
}

Expr Group::subs(const Expr::subst_map& map) const
{
  return Group(expr.subs(map)).clone();
}

Expr Group::integrate_internal(const Symbol& s) const
{
  return Group(expr.integrate(s)).clone();
}

Expr Group::simplify() const
{
  const Expr simplified = expr.simplify();

  if (!is_a<Rational>(expr.internal()) && !is_a<Float>(expr.internal()))
  {
    return Group(simplified).clone();
  }
  else
  {
    return simplified;
  }
}

Float Group::eval(const Expr::subst_map& map) const
{
  return expr.eval(map);
}

Expr Group::getExpr() const
{
  return expr;
}

bool Group::operator==(const Group& g) const
{
  return expr == g.expr;
}

std::size_t Group::untypedHash() const
{
  std::size_t result = 0x2f33cf1c;
  boost::hash_combine(result, expr.hashValue());
  return result;
}

void Group::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr.accept(v);
}

}

}
