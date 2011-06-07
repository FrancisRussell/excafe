#ifndef SIMPLE_CFD_NUMERIC_EXPRESSION_VISITOR_HPP
#define SIMPLE_CFD_NUMERIC_EXPRESSION_VISITOR_HPP

#include <cln/cln.h>

namespace cfd
{

template<typename V>
class NumericExpressionVisitor
{
public:
  typedef V         variable_t;
  typedef cln::cl_R value_t;

  virtual void visitConstant(const value_t& s) = 0;
  virtual void visitVariable(const variable_t& var) = 0;
  virtual void visitExponent(const int exponent) = 0;
  virtual void visitAbsoluteValue() = 0;
  virtual void postSummation(const std::size_t nops) = 0;
  virtual void postProduct(const std::size_t nops) = 0;

  virtual ~NumericExpressionVisitor() {}
};

}

#endif
