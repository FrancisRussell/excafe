#include <simple_cfd/symbolic/group.hpp>
#include <simple_cfd/symbolic/float.hpp>
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
  o << expr;
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

Expr Group::integrate(const Symbol& s) const
{
  return Group(expr.integrate(s)).clone();
}

Expr Group::simplify() const
{
  return Group(expr.simplify()).clone();
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
  return expr.hashValue();
}

void Group::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr.accept(v);
}

}

}
