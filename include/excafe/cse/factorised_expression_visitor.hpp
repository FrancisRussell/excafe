#ifndef EXCAFE_CSE_FACTORISED_EXPRESSION_VISITOR_HPP
#define EXCAFE_CSE_FACTORISED_EXPRESSION_VISITOR_HPP

#include "polynomial_index.hpp"
#include <excafe/numeric/expression_visitor.hpp>

namespace excafe
{

namespace cse
{

template<typename V>
class FactorisedExpressionVisitor : public NumericExpressionVisitor<V>
{
public:
  virtual void visitOriginalTerm(const unsigned index) = 0;
  virtual void visitFactorisedTerm(const PolynomialIndex& i) = 0;
  virtual void postOriginalTerm(const unsigned index) = 0;
  virtual void postFactorisedTerm(const PolynomialIndex& i) = 0;
};

}

}

#endif
