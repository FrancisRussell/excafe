#ifndef SIMPLE_CFD_SYMBOLIC_EXPR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/operators.hpp>
#include <map>
#include "symbolic_fwd.hpp"
#include <simple_cfd/numeric/expression.hpp>
#include <simple_cfd/numeric/expression_visitor.hpp>

namespace cfd
{

namespace symbolic
{

class Expr : public NumericExpression<Symbol>, 
                    boost::arithmetic<Expr>
{
public:
  typedef boost::shared_ptr<const Basic> ref_t;

private:
  ref_t expr;

public:
  typedef std::map<Symbol, Expr> subst_map;

  explicit Expr(Basic* e);
  explicit Expr(ref_t& e);
  Expr();
  Expr(const double s);
  Expr& operator=(const Expr& e);
  bool operator==(const Expr& e) const;
  bool operator!=(const Expr& e) const;
  Expr& operator+=(const Expr& e);
  Expr& operator-=(const Expr& e);
  Expr& operator/=(const Expr& e);
  Expr& operator*=(const Expr& e);
  Expr operator-() const;
  bool has(const Expr& e) const;
  void write(std::ostream& o) const;
  std::size_t hashValue() const;
  Expr derivative(const Symbol& s) const;
  Expr simplify() const;
  Expr integrate(const Symbol& s) const;
  Expr integrate(const Symbol& s, const Float& a, const Float& b) const;
  const Basic& internal() const;
  Expr subs(const subst_map& map) const;
  void accept(Visitor& v) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  void traverse(Visitor& v) const;
  void swap(Expr& e);
  Expr expand() const;
  Expr eval() const;
};

std::size_t hash_value(const Expr& e);
Expr pow(const Expr& e, int power);
std::ostream& operator<<(std::ostream& o, const Expr& e);

}

}


#endif
