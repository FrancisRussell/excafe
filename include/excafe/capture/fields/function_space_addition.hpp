#ifndef EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_ADDITION_HPP
#define EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_ADDITION_HPP

#include "function_space_binary_operator.hpp"
#include <boost/shared_ptr.hpp>

namespace excafe
{

namespace detail
{

class FunctionSpaceAddition : public FunctionSpaceBinaryOperator
{
public:
  FunctionSpaceAddition(const FunctionSpaceExpr::expr_ptr& _left,
                        const FunctionSpaceExpr::expr_ptr& _right) :
    FunctionSpaceBinaryOperator(_left, _right)
  {
  }

  void accept(FunctionSpaceVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
