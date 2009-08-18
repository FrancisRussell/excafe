#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_EMPTY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_EMPTY_HPP

#include "function_space_expr.hpp"
#include "function_space_visitor.hpp"

namespace cfd
{

namespace detail
{

class FunctionSpaceEmpty : public FunctionSpaceExpr
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
