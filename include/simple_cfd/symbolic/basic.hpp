#ifndef SIMPLE_CFD_SYMBOLIC_BASIC_HPP
#define SIMPLE_CFD_SYMBOLIC_BASIC_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <boost/shared_ptr.hpp>
#include "symbolic_fwd.hpp"
#include "expr.hpp"
#include <simple_cfd/util/type_info.hpp>

namespace cfd
{

namespace symbolic
{

class Basic
{
public:
  virtual std::size_t nops() const = 0;
  virtual void write(std::ostream& o) const = 0;
  virtual Expr clone() const = 0;
  virtual Expr derivative(const Symbol& s) const = 0;
  virtual Expr integrate(const Symbol& s) const = 0;
  virtual Expr simplify() const = 0;
  virtual bool isNumber() const = 0;
  virtual std::size_t hashValue() const = 0;
  virtual bool operator==(const Basic& b) const = 0;
  virtual bool operator<(const Basic& b) const = 0;
  virtual operator Expr() const = 0;
  virtual Expr subs(const Expr::subst_map& map) const = 0;
  virtual bool has(const Expr& e) const = 0;
  virtual void accept(Visitor& v) const = 0;
  virtual ~Basic() {}
};

std::ostream& operator<<(std::ostream& o, const Basic& b);

}

}

#endif
