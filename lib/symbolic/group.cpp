#include <excafe/symbolic/group.hpp>
#include <excafe/symbolic/float.hpp>
#include <excafe/symbolic/rational.hpp>
#include <excafe/symbolic/symbol.hpp>
#include <excafe/symbolic/expr.hpp>
#include <excafe/symbolic/extracted_expressions.hpp>
#include <excafe/util/hash.hpp>
#include <set>
#include <ostream>
#include <cstddef>

namespace excafe
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

Expr Group::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  return Group(expr.derivative(s, cache)).clone();
}

bool Group::depends(const std::set<Symbol>& symbols) const
{
  return expr.depends(symbols);
}

Expr Group::subs(const Expr::subst_map& map, const unsigned flags) const
{
  return Group(expr.subs(map, flags)).clone();
}

Expr Group::integrate(const Symbol& s, const unsigned flags) const
{
  return Group(expr.integrate(s, flags)).clone();
}

Expr Group::simplify() const
{
  const Expr simplified = expr.simplify();

  if (!is_exactly_a<Rational>(simplified) && !is_exactly_a<Float>(simplified))
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
  excafe::util::hash_accum(result, expr.hashValue());
  return result;
}

void Group::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr.accept(v);
}

Expr Group::extractPolynomials(ExtractedExpressions& extracted) const
{
  const Symbol subExpr = extracted.addRepresentable(expr.extractPolynomials(extracted));
  return extracted.addUnrepresentable(Group(subExpr));
}

}

}
