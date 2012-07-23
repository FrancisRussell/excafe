#include <excafe/symbolic/abs.hpp>
#include <excafe/symbolic/float.hpp>
#include <excafe/symbolic/rational.hpp>
#include <excafe/symbolic/expr.hpp>
#include <excafe/symbolic/product.hpp>
#include <excafe/symbolic/symbol.hpp>
#include <excafe/symbolic/extracted_expressions.hpp>
#include <excafe/util/hash.hpp>
#include <excafe/exception.hpp>
#include <set>
#include <ostream>
#include <cstddef>

namespace excafe
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
  excafe::util::hash_accum(result, expr.hashValue());
  return result;
}

void Abs::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr.accept(v);
  v.visitAbsoluteValue();
}

bool Abs::isPolynomial() const
{
  return false;
}

Expr Abs::extractPolynomials(ExtractedExpressions& extracted) const
{
  const Symbol subExpr = extracted.addRepresentable(expr.extractPolynomials(extracted));
  return extracted.addUnrepresentable(Abs(subExpr));
}

}

}
