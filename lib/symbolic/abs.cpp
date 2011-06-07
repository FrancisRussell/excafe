#include <simple_cfd/symbolic/abs.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/util/hash.hpp>
#include <simple_cfd/exception.hpp>
#include <set>
#include <ostream>
#include <cstddef>

namespace cfd
{

namespace symbolic
{

Abs::Abs(const Expr& e) : expr(e)
{
}

std::size_t Abs::nops() const
{
  return 1;
}

void Abs::write(std::ostream& o) const
{
  o << "abs(" << expr << ")";
}

Expr Abs::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  return Product::mul(sign(), expr.derivative(s, cache)).clone();
}

bool Abs::depends(const std::set<Symbol>& symbols) const
{
  return expr.depends(symbols);
}

Expr Abs::subs(const Expr::subst_map& map, const unsigned flags) const
{
  return Abs(expr.subs(map, flags)).clone();
}

Expr Abs::integrate(const Symbol& s, const unsigned flags) const
{
  CFD_EXCEPTION("integral of abs() unimplemented.");
}

Expr Abs::simplify() const
{
  const Expr simplified = expr.simplify();

  if (!is_exactly_a<Rational>(simplified))
  {
    return Abs(simplified).clone();
  }
  else
  {
    return convert_to<Rational>(simplified).abs();
  }
}

Float Abs::eval(const Expr::subst_map& map) const
{
  return expr.eval(map).abs();
}

Expr Abs::getExpr() const
{
  return expr;
}

Expr Abs::sign() const
{
  return Product::div(expr, *this);
}

bool Abs::operator==(const Abs& a) const
{
  return expr == a.expr;
}

std::size_t Abs::untypedHash() const
{
  std::size_t result = 0x2371d80b;
  cfd::util::hash_accum(result, expr.hashValue());
  return result;
}

void Abs::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr.accept(v);
  v.visitAbsoluteValue();
}

}

}
