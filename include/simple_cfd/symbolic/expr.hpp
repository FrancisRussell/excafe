#ifndef SIMPLE_CFD_SYMBOLIC_EXPR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPR_HPP

#include <boost/shared_ptr.hpp>
#include "symbolic_fwd.hpp"

namespace cfd
{

namespace symbolic
{

class Expr
{
public:
  typedef boost::shared_ptr<const Basic> ref_t;

private:
  ref_t expr;

public:
  explicit Expr(Basic* e);
  explicit Expr(ref_t& e);

  Expr& operator=(const Expr& e);
  
  int getTypeID() const;
};

Expr operator+(const Expr& a, const Expr& b);

}

}


#endif
