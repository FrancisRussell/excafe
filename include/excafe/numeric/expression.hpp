#ifndef EXCAFE_NUMERIC_EXPRESSION_HPP
#define EXCAFE_NUMERIC_EXPRESSION_HPP

#include "numeric_fwd.hpp"

namespace excafe
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
