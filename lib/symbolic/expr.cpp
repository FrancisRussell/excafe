#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/basic.hpp>
#include <simple_cfd/symbolic/number.hpp>
#include <simple_cfd/symbolic/sum.hpp>
#include <simple_cfd/symbolic/product.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expand_visitor.hpp>
#include <ostream>

namespace cfd
{

namespace symbolic
{

Expr::Expr(Basic* e) : expr(e)
{
}

Expr::Expr(ref_t& e) : expr(e)
{
}

Expr::Expr(const double s) : expr(new Number(s))
{
}

Expr::Expr() : expr(new Number(0))
{
}

Expr operator+(const Expr& a, const Expr& b)
{
  return Sum(a, b).simplify();
}

Expr operator-(const Expr& a, const Expr& b)
{
  return Sum::sub(a, b).simplify();
}

Expr operator*(const Expr& a, const Expr& b)
{
  return Product(a, b).simplify();
}

Expr operator/(const Expr& a, const Expr& b)
{
  return Product::div(a, b).simplify();
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

bool Expr::operator<(const Expr& e) const
{
  const bool result = *expr < *e.expr;
  return result;
}

bool Expr::operator==(const Expr& e) const
{
  const bool result = *expr == *e.expr;
  return result;
}

bool Expr::operator!=(const Expr& e) const
{
  return !(*this == e);
}

Expr Expr::operator-() const
{
  return Sum::integer_multiple(*this, -1);
}

bool Expr::has(const Expr& e) const
{
  return expr->has(e);
}

std::size_t Expr::hashValue() const
{
  return expr->hashValue();
}

Expr Expr::simplify() const
{
  return expr->simplify();
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

Expr Expr::derivative(const Symbol& s) const
{
  return expr->derivative(s).simplify();
}

Expr Expr::integrate(const Symbol& s) const
{
  return expr->integrate(s).simplify();
}

Expr Expr::integrate(const Symbol& s, const Number& a, const Number& b) const
{
  const Expr integrated = expr->integrate(s);

  subst_map aMap, bMap;
  aMap[s] = a;
  bMap[s] = b;

  return (integrated.subs(bMap) - integrated.subs(aMap)).simplify();
}

void Expr::accept(Visitor& v) const
{
  return expr->accept(v);
}

Expr Expr::expand() const
{
  ExpandVisitor v;
  expr->accept(v);
  return v.getResult().simplify();
}

void Expr::accept(NumericExpressionVisitor<Symbol>& v) const
{
  expr->accept(v);
}

Expr Expr::subs(const subst_map& map) const
{
  return expr->subs(map).simplify();
}

Expr pow(const Expr& e, const int power)
{
  return Product::pow(e, power).simplify();
}

}

}
