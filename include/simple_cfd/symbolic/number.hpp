#ifndef SIMPLE_CFD_SYMBOLIC_NUMBER_HPP
#define SIMPLE_CFD_SYMBOLIC_NUMBER_HPP

#include <string>
#include <cstddef>
#include <ostream>
#include "symbolic_fwd.hpp"
#include "abstract_basic.hpp"

namespace cfd
{

namespace symbolic
{

class Number : public AbstractBasic<Number>
{
private:
  double value;

public:
  Number(const double _value);
  virtual std::size_t nops() const;
  virtual void write(std::ostream& o) const;
  virtual Expr derivative(const Symbol& s) const;
  virtual bool isNumber() const;
  bool operator==(const Number& n) const;
  bool operator<(const Number& n) const;
  std::size_t untypedHash() const;
};

}

}

#endif
