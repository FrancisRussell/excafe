#ifndef SIMPLE_CFD_SYMBOLIC_PAIR_SET_HPP
#define SIMPLE_CFD_SYMBOLIC_PAIR_SET_HPP

#include <cstddef>
#include <ostream>
#include <utility>
#include <vector>
#include "basic.hpp"

namespace cfd
{

namespace symbolic
{

class Sum : public Basic
{
private:
  std::set<Basic::ref_t> terms;

public:
  virtual std::size_t nops() const
  {
    return terms.size();
  }

  virtual void write(std::ostream& o) const
  {
  }

  virtual Expr clone() const = 0;
  virtual Expr derivative(const Symbol& s) const = 0;
  virtual bool isNumber() const = 0;
  virtual ~Basic();
};

}

}

#endif
