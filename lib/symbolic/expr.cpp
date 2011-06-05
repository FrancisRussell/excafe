#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/basic.hpp>
#include <simple_cfd/symbolic/float.hpp>
#include <simple_cfd/symbolic/rational.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expand_visitor.hpp>
#include <simple_cfd/symbolic/make_expr_from.hpp>
#include <simple_cfd/symbolic/flags.hpp>
#include <ostream>
#include <cassert>
#include <boost/optional.hpp>

namespace cfd
{

namespace symbolic
{

Expr Expr::initial = make_expr_from(Rational(0));

Expr::Expr(const ref_t& e) : expr(e)
{
  // e should already been marked as heap allocated. We cannot set it since e is const.
}

Expr::Expr(const double s) : expr(make_expr_from(Float(s)).expr)
{
}

Expr::Expr(const long s) : expr(make_expr_from(Rational(s)).expr)
{
}

// The default contructor uses a shared initial value to avoid invoking
// malloc on each construction.
Expr::Expr() : expr(initial.expr)
{
}

Expr& Expr::operator+=(const Expr& e)
{
  Expr newExpr = Sum::add(*this, e).simplify();
  std::swap(newExpr, *this);
  return *this;
}

Expr& Expr::operator-=(const Expr& e)
{
  Expr newExpr = Sum::sub(*this, e).simplify();
  std::swap(newExpr, *this);
  return *this;
}

Expr& Expr::operator*=(const Expr& e)
{
  Expr newExpr = Product::mul(*this, e).simplify();
  std::swap(newExpr, *this);
  return *this;
}

Expr& Expr::operator/=(const Expr& e)
{
  Expr newExpr = Product::div(*this, e).simplify();
  std::swap(newExpr, *this);
  return *this;
}

Expr& Expr::operator=(const Expr& e)
{
  expr = e.expr;
  return *this;
}

void Expr::write(std::ostream& o) const
{
  o << *expr;
}

bool Expr::operator==(const Expr& e) const
{
  if (expr.get() == e.expr.get())
    return true;
  else
    return *expr == *e.expr;
}

Expr Expr::operator-() const
{
  return Sum::rational_multiple(*this, Rational(-1)).simplify();
}

bool Expr::depends(const std::set<Symbol>& symbols) const
{
  return expr->depends(symbols);
}

bool Expr::depends(const Symbol& s) const
{
  std::set<Symbol> symbols;
  symbols.insert(s);
  return depends(symbols);
}

std::size_t Expr::hashValue() const
{
  return expr->hashValue();
}

Expr Expr::simplify() const
{
  return expr->simplify();
}

std::size_t hash_value(const Basic& b)
{
  return b.hashValue();
}

std::size_t hash_value(const Expr& e)
{
  return e.hashValue();
}

const Basic& Expr::internal() const
{
  return *expr;
}

std::ostream& operator<<(std::ostream& o, const Expr& e)
{
  e.write(o);
  return o;
}

Expr Expr::derivative(const Symbol& s, Expr::optional_derivative_cache cache) const
{
  if (cache)
    return cache->derivative(*this, s);
  else
    return expr->derivative(s, cache).simplify();
}

Expr Expr::integrate(const Symbol& s, const unsigned flags) const
{
  const Expr integrated = expr->integrate(s, flags);
  return integrated.simplify();
}

Expr Expr::integrate(const region_t& region, const unsigned flags) const
{
  const Expr integrated = expr->integrate(region, flags);
  return integrated.simplify();
}

Expr Expr::integrate(const Symbol& s, const Rational& a, const Rational& b, const unsigned flags) const
{
  region_t region;
  region.setInterval(s, a, b);
  return this->integrate(region, flags);
}

void Expr::accept(Visitor& v) const
{
  return expr->accept(v);
}

Expr Expr::expand(const bool distribute) const
{
  ExpandVisitor v;
  expr->accept(v);
  return v.getResult(distribute).simplify();
}

void Expr::traverse(Visitor& v) const
{
  expr->accept(v);
}

void Expr::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr->accept(v);
}

Expr Expr::subs(const subst_map& subs, const unsigned flags) const
{
  Expr result;

  if ((flags & Flags::DO_NOT_SIMPLIFY_SUBST_MAP) == 0)
  {
    subst_map simplifiedSubs;
    BOOST_FOREACH(const subst_map::value_type& sub, subs)
      simplifiedSubs.insert(subst_map::value_type(sub.first, sub.second.simplify()));

    return expr->subs(simplifiedSubs, flags | Flags::DO_NOT_SIMPLIFY_SUBST_MAP).simplify();
  }
  else
  {
    return expr->subs(subs, flags).simplify();
  }
}

Expr pow(const Expr& e, const int power)
{
  return Product::pow(e, power).simplify();
}

Float Expr::eval(const subst_map& map) const
{
  return expr->eval(map);
}

}

}
