#ifndef SIMPLE_CFD_NUMERIC_EXPRESSION_HPP
#define SIMPLE_CFD_NUMERIC_EXPRESSION_HPP

#include "numeric_fwd.hpp"

namespace cfd
{

template<typename V>
class NumericExpression
{
public:
  virtual void accept(NumericExpressionVisitor<V>& v) const = 0;
  virtual ~NumericExpression() {}
};

}

#endif
