#ifndef SIMPLE_CFD_SYMBOLIC_EXPR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/operators.hpp>
#include <map>
#include <set>
#include <cassert>
#include "symbolic_fwd.hpp"
#include <simple_cfd/numeric/expression.hpp>
#include <simple_cfd/numeric/expression_visitor.hpp>

namespace cfd
{

namespace symbolic
{

class Expr : public NumericExpression<Symbol>, 
                    boost::arithmetic<Expr>,
                    boost::equality_comparable<Expr>
{
public:
  typedef boost::shared_ptr<const Basic> ref_t;

private:
  template<typename T> friend bool is_a(const Expr&);
  template<typename T> friend const T& convert_to(const Expr&);
  static Expr initial;
  ref_t expr;

public:
  typedef std::map<Symbol, Expr> subst_map;

  explicit Expr(Basic* e);
  explicit Expr(const ref_t& e);
  Expr();
  Expr(const double s);
  Expr& operator=(const Expr& e);
  bool operator==(const Expr& e) const;
  Expr& operator+=(const Expr& e);
  Expr& operator-=(const Expr& e);
  Expr& operator/=(const Expr& e);
  Expr& operator*=(const Expr& e);
  Expr operator-() const;
  bool depends(const Symbol& s) const;
  bool depends(const std::set<Symbol>& symbols) const;
  void write(std::ostream& o) const;
  std::size_t hashValue() const;
  Expr derivative(const Symbol& s) const;
  Expr simplify() const;
  Expr integrate(const Symbol& s) const;
  Expr integrate_internal(const Symbol& s) const;
  Expr integrate(const Symbol& s, const Float& a, const Float& b) const;
  const Basic& internal() const;
  Expr subs(const subst_map& map) const;
  void accept(Visitor& v) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  void traverse(Visitor& v) const;
  void swap(Expr& e);
  Expr expand() const;
  Float eval(const subst_map& map) const;
};

template<typename T>
inline bool is_a(const Basic& e)
{
  return dynamic_cast<const T*>(&e) != NULL;
}

template<typename T>
inline bool is_exactly_a(const Basic& e)
{
  return typeid(e) == typeid(T);
}

template<typename T>
inline const T& convert_to(const Basic& e)
{
  return static_cast<const T&>(e);
}

template<typename T>
inline bool is_a(const Expr& e)
{
  return is_a<T>(e.internal());
}

template<typename T>
inline bool is_exactly_a(const Expr& e)
{
  return is_exactly_a<T>(e.internal());
}

template<typename T>
inline const T& convert_to(const Expr& e)
{
  return convert_to<T>(e.internal());
}

std::size_t hash_value(const Expr& e);
std::size_t hash_value(const Basic& b);

Expr pow(const Expr& e, int power);
std::ostream& operator<<(std::ostream& o, const Expr& e);

}

}


#endif
