#ifndef SIMPLE_CFD_NUMERIC_EXPRESSION_VISITOR_HPP
#define SIMPLE_CFD_NUMERIC_EXPRESSION_VISITOR_HPP

#include <cln/cln.h>
#include <simple_cfd/mp/integer.hpp>

namespace cfd
{

template<typename V>
class NumericExpressionVisitor
{
public:
  typedef V          variable_t;
  typedef cln::cl_F   float_t;
  typedef mp::Integer integer_t;

  virtual void visitConstant(const float_t& s) = 0;
  virtual void visitConstant(const integer_t& s) = 0;
  virtual void visitVariable(const variable_t& var) = 0;
  virtual void visitExponent(const int exponent) = 0;
  virtual void visitAbsoluteValue() = 0;
  virtual void postSummation(const std::size_t nops) = 0;
  virtual void postProduct(const std::size_t nops) = 0;

  virtual ~NumericExpressionVisitor() {}
};

}

#endif
