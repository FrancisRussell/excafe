#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_BINARY_OPERATOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_BINARY_OPERATOR_HPP

#include <boost/shared_ptr.hpp>
#include "function_space_expr.hpp"

namespace cfd
{

namespace detail
{

class FunctionSpaceBinaryOperator : public FunctionSpaceExpr
{
private:
  FunctionSpaceExpr::expr_ptr left;
  FunctionSpaceExpr::expr_ptr right;

public:
  FunctionSpaceBinaryOperator(const FunctionSpaceExpr::expr_ptr& _left,
                              const FunctionSpaceExpr::expr_ptr& _right) :
    left(_left), right(_right)
  {
  }

  FunctionSpaceExpr& getLeft()
  {
    return *left;
  }

  FunctionSpaceExpr& getRight()
  {
    return *right;
  }
};

}

}

#endif
