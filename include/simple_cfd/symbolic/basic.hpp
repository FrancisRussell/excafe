#ifndef SIMPLE_CFD_SYMBOLIC_BASIC_HPP
#define SIMPLE_CFD_SYMBOLIC_BASIC_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include <boost/shared_ptr.hpp>
#include <boost/mpl/void.hpp>
#include "symbolic_fwd.hpp"
#include "ertti.hpp"
#include "expr.hpp"

namespace cfd
{

namespace symbolic
{

class Basic : public ERTTI<Basic, Basic>
{
public:
  virtual std::size_t nops() const = 0;
  virtual void write(std::ostream& o) const = 0;
  virtual Expr clone() const = 0;
  virtual Expr derivative(const Symbol& s) const = 0;
  virtual bool isNumber() const = 0;
  virtual std::size_t hashValue() const = 0;
  virtual ~Basic() {}

  operator Expr() const
  {
    return Expr(clone());
  }
};

}

}

#endif
