#ifndef EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_UNDEFINED_HPP
#define EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_UNDEFINED_HPP

#include "function_space_expr.hpp"
#include "function_space_visitor.hpp"

namespace excafe
{

namespace detail
{

class FunctionSpaceUndefined : public FunctionSpaceExpr
{
public:
  void accept(FunctionSpaceVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
