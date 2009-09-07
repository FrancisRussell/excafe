#include <simple_cfd/capture/fields/function_space.hpp>
#include <simple_cfd/capture/fields/function_space_undefined.hpp>
#include <simple_cfd/capture/fields/function_space_addition.hpp>


namespace cfd
{

FunctionSpace::FunctionSpace() : expr(new detail::FunctionSpaceUndefined())
{
}

FunctionSpace::FunctionSpace(detail::FunctionSpaceExpr* const _expr) : expr(_expr)
{
}

FunctionSpace::FunctionSpace& FunctionSpace::operator+=(const FunctionSpace& f)
{
  expr = expr_ptr(new detail::FunctionSpaceAddition(expr, f.expr));
  return *this;
}

FunctionSpace::expr_ptr FunctionSpace::getExpr() const
{
    return expr;
}

}
