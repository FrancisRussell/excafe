#ifndef SIMPLE_CFD_SYMBOLIC_BASIC_HPP
#define SIMPLE_CFD_SYMBOLIC_BASIC_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <boost/shared_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>
#include "symbolic_fwd.hpp"
#include "expr.hpp"
#include <simple_cfd/util/type_info.hpp>
#include <simple_cfd/numeric/expression_visitor.hpp>

namespace cfd
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
  virtual Expr derivative(const Symbol& s) const = 0;
  virtual Expr integrate_internal(const Symbol& s) const = 0;
  virtual Expr simplify() const = 0;
  virtual Expr extractMultiplier(Rational& coeff) const = 0;
  virtual Float eval(const Expr::subst_map& map) const = 0;
  virtual std::size_t hashValue() const = 0;
  virtual bool operator==(const Basic& b) const = 0;
  virtual operator Expr() const = 0;
  virtual Expr subs(const Expr::subst_map& map) const = 0;
  virtual bool has(const Expr& e) const = 0;
  virtual void accept(Visitor& v) const = 0;
  virtual void accept(NumericExpressionVisitor<Symbol>& v) const = 0;
  virtual ~Basic() {}
};

std::ostream& operator<<(std::ostream& o, const Basic& b);

}

}

#endif
