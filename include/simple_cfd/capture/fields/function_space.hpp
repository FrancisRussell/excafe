#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_HPP

#include "fields_fwd.hpp"
#include "function_space_expr.hpp"
#include "function_space_empty.hpp"
#include "function_space_addition.hpp"
#include <boost/operators.hpp>

namespace cfd
{

class FunctionSpace : boost::addable<FunctionSpace>
{
private:
  typedef detail::FunctionSpaceExpr::expr_ptr expr_ptr;
  expr_ptr expr;

public:
  FunctionSpace() : expr(new detail::FunctionSpaceEmpty())
  {
  }

  FunctionSpace(detail::FunctionSpaceExpr* const _expr) : expr(_expr)
  {
  }
  
  FunctionSpace& operator+=(const FunctionSpace& f)
  {
    expr = expr_ptr(new detail::FunctionSpaceAddition(expr, f.expr)); 
    return *this;
  }
};

}

#endif
