#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP

#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>

namespace cfd
{

namespace detail
{

namespace
{

class FormTemporalIndexFinder : public FieldVisitor
{
private:
  TemporalIndexSet indices;

public:
  TemporalIndexSet getIndices() const
  {
    return indices;
  }

  virtual void enter(FieldAddition& addition) {}
  virtual void exit(FieldAddition& addition) {}

  virtual void enter(FieldInnerProduct& inner) {}
  virtual void exit(FieldInnerProduct& inner) {}

  virtual void enter(FieldOuterProduct& outer) {}
  virtual void exit(FieldOuterProduct& outer) {}

  virtual void enter(FieldColonProduct& colon) {}
  virtual void exit(FieldColonProduct& colon) {}

  virtual void enter(FieldGradient& gradient) {}
  virtual void exit(FieldGradient& gradient) {}

  virtual void enter(FieldDivergence& divergence) {}
  virtual void exit(FieldDivergence& divergence) {}

  // Terminals
  virtual void visit(FacetNormal& normal) {}
  virtual void visit(FieldBasis& basis) {}
  virtual void visit(FieldDiscreteReference& field)
  {
    indices += field.getDiscreteField().getExpr()->getTemporalIndices();
  }

  virtual void visit(FieldScalar& s)
  {
    indices += s.getValue().getExpr()->getTemporalIndices();
  }
};

}

class OperatorAssembly : public OperatorExpr
{
private:
  const FunctionSpaceExpr::expr_ptr trialSpace;
  const FunctionSpaceExpr::expr_ptr testSpace;
  const forms::BilinearFormIntegralSum sum;

  template<typename ForwardIterator>
  TemporalIndexSet getTemporalIndices(const ForwardIterator begin, const ForwardIterator end) const
  {
    FormTemporalIndexFinder indexFinder;

    for(ForwardIterator i(begin); i!=end; ++i)
    {
      i->getTrialField()->accept(indexFinder);
      i->getTestField()->accept(indexFinder);
    }
    
    return indexFinder.getIndices();
  }

public:
  OperatorAssembly(const FunctionSpaceExpr::expr_ptr& _trialSpace, const FunctionSpaceExpr::expr_ptr& _testSpace,
                   const forms::BilinearFormIntegralSum& _sum) : 
    trialSpace(_trialSpace), testSpace(_testSpace), sum(_sum)
  {
  }

  forms::BilinearFormIntegralSum getBilinearFormIntegralSum() const
  {
    return sum;
  }

  void accept(DiscreteExprVisitor& v)
  {
    v.visit(*this);
  }

  virtual FunctionSpaceExpr::expr_ptr getTrialSpace() const
  {
    return trialSpace;
  }

  virtual FunctionSpaceExpr::expr_ptr getTestSpace() const
  {
    return testSpace;
  }

  virtual TemporalIndexSet getTemporalIndices() const
  {
    TemporalIndexSet indices;
    indices += getTemporalIndices(sum.begin_dx(), sum.end_dx());
    indices += getTemporalIndices(sum.begin_ds(), sum.end_ds());
    indices += getTemporalIndices(sum.begin_dS(), sum.end_dS());

    return indices;
  }
};

}

}

#endif
