#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_HPP

#include "fields_fwd.hpp"
#include "function_space_expr.hpp"
#include <boost/operators.hpp>

namespace cfd
{

class FunctionSpace : boost::addable<FunctionSpace>
{
public:
  typedef detail::FunctionSpaceExpr::expr_ptr expr_ptr;

private:
  expr_ptr expr;

public:
  FunctionSpace();
  FunctionSpace(detail::FunctionSpaceExpr* const _expr);
  FunctionSpace(const detail::FunctionSpaceExpr::expr_ptr _expr);
  FunctionSpace& operator+=(const FunctionSpace& f);
  expr_ptr getExpr() const;
};

}

#endif
