#ifndef SIMPLE_CFD_CSE_FACTORISED_EXPRESSION_VISITOR_HPP
#define SIMPLE_CFD_CSE_FACTORISED_EXPRESSION_VISITOR_HPP

#include "polynomial_index.hpp"
#include <simple_cfd/numeric/expression_visitor.hpp>

namespace cfd
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
