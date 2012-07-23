#ifndef EXCAFE_SYMBOLIC_EXPR_HPP
#define EXCAFE_SYMBOLIC_EXPR_HPP

#include <boost/shared_ptr.hpp>
#include <boost/operators.hpp>
#include <boost/optional.hpp>
#include <map>
#include <set>
#include <cassert>
#include "symbolic_fwd.hpp"
#include "derivative_cache.hpp"
#include <excafe/numeric/expression.hpp>
#include <excafe/numeric/expression_visitor.hpp>
#include <excafe/numeric/orthotope.hpp>
#include <cln/float.h>
#include <cln/rational.h>

namespace excafe
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
  static Expr initial;
  ref_t expr;

public:
  static const bool supports_abs = true;

  typedef Symbol variable_t;
  typedef std::map<Symbol, Expr> subst_map;
  typedef numeric::Orthotope<Symbol, Rational> region_t;
  typedef DerivativeCache derivative_cache;
  typedef boost::optional<derivative_cache&> optional_derivative_cache;

  explicit Expr(const ref_t& e);
  Expr();
  Expr(const double s);
  Expr(const long s);
  Expr(const int s);
  Expr(const cln::cl_RA& s);
  Expr(const cln::cl_F& s);
  Expr& operator=(const Expr& e);
  bool operator==(const Expr& e) const;
  Expr& operator+=(const Expr& e);
  Expr& operator-=(const Expr& e);
  Expr& operator/=(const Expr& e);
  Expr& operator*=(const Expr& e);
  Expr operator-() const;
  bool depends(const Symbol& s) const;
  bool depends(const std::set<Symbol>& symbols) const;
  std::set<Symbol> getSymbols() const;
  void write(std::ostream& o) const;
  std::size_t hashValue() const;
  Expr derivative(const Symbol& s, 
                  optional_derivative_cache cache = optional_derivative_cache()) const;
  Expr simplify() const;
  Expr extractMultiplier(Rational& r) const;
  Expr integrate(const Symbol& s, unsigned flags = 0) const;
  Expr integrate(const Symbol& s, const Rational& a, const Rational& b, unsigned flags = 0) const;
  Expr integrate(const region_t& region, unsigned flags = 0) const;
  const Basic& internal() const;
  Expr subs(const subst_map& map, unsigned flags = 0) const;
  void accept(Visitor& v) const;
  void accept(NumericExpressionVisitor<Symbol>& v) const;
  void traverse(Visitor& v) const;
  void swap(Expr& e);
  Expr expand(bool distribute = true) const;
  Expr extractPolynomials(ExtractedExpressions& extracted) const;
  bool isPolynomial() const;
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
Expr abs(const Expr& e);
std::ostream& operator<<(std::ostream& o, const Expr& e);

}

}


#endif
