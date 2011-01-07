#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/basic.hpp>

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

Expr operator+(const Expr& a, const Expr& b)
{
  return a;
}

Expr& Expr::operator=(const Expr& e)
{
  expr = e.expr;
  return *this;
}

int Expr::getTypeID() const
{
  return expr->getTypeID();
}

}

}
