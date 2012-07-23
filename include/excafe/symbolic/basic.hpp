#ifndef EXCAFE_SYMBOLIC_BASIC_HPP
#define EXCAFE_SYMBOLIC_BASIC_HPP

#include <string>
#include <cstddef>
#include <set>
#include <ostream>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "symbolic_fwd.hpp"
#include "expr.hpp"
#include <excafe/util/type_info.hpp>
#include <excafe/numeric/expression_visitor.hpp>

namespace excafe
{

namespace symbolic
{

class Basic : public boost::enable_shared_from_this<Basic>
{
public:
  virtual void markHeapAllocated() = 0;
  virtual std::size_t nops() const = 0;
  virtual void write(std::ostream& o) const = 0;
  virtual Expr clone() const = 0;
  virtual Expr derivative(const Symbol& s, Expr::optional_derivative_cache cache) const = 0;
  virtual Expr integrate(const Symbol& s, unsigned flags) const = 0;
  virtual Expr integrate(const Expr::region_t& region, unsigned flags) const = 0;
  virtual Expr simplify() const = 0;
  virtual Expr extractMultiplier(Rational& coeff) const = 0;
  virtual Float eval(const Expr::subst_map& map) const = 0;
  virtual std::size_t typeHash() const = 0;
  virtual std::size_t hashValue() const = 0;
  virtual bool operator==(const Basic& b) const = 0;
  virtual operator Expr() const = 0;
  virtual Expr subs(const Expr::subst_map& map, unsigned flags) const = 0;
  virtual bool depends(const std::set<Symbol>& e) const = 0;
  virtual void accept(Visitor& v) const = 0;
  virtual void accept(NumericExpressionVisitor<Symbol>& v) const = 0;
  virtual bool isPolynomial() const = 0;
  virtual Expr extractPolynomials(ExtractedExpressions& extracted) const = 0;
  virtual ~Basic() {}
};

std::ostream& operator<<(std::ostream& o, const Basic& b);

}

}

#endif
