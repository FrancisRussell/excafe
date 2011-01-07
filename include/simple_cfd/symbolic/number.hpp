#ifndef SIMPLE_CFD_SYMBOLIC_NUMBER_HPP
#define SIMPLE_CFD_SYMBOLIC_NUMBER_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include "symbolic_fwd.hpp"
#include "basic.hpp"

namespace cfd
{

namespace symbolic
{

DECLARE_SYMBOLIC_NODE(Number, Basic)
{
private:
  double value;

public:
  Number(const double _value);

  virtual Expr clone() const;

  virtual std::size_t nops() const;

  virtual void write(std::ostream& o) const;

  virtual Expr derivative(const Symbol& s) const;

  virtual bool isNumber() const;

  std::size_t hashValue() const;
};

}

}

#endif
