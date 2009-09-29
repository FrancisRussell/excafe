#ifndef SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP

#include <memory>
#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include <simple_cfd/capture/forms/bilinear_form_integral_sum.hpp>
#include <simple_cfd/capture/indices/propagation_rules.hpp>

namespace cfd
{

namespace detail
{

class FormPropagationRulesGetter : public FieldVisitor
{
private:
  OperatorAssembly& parent;
  PropagationRules rules;

public:
  FormPropagationRulesGetter(OperatorAssembly& _parent) : parent(_parent)
  {
  }
  
  PropagationRules getRules() const;

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

  // References to discrete expressions
  virtual void visit(FieldDiscreteReference& field);
  virtual void visit(FieldScalar& s);
};

class OperatorAssembly : public OperatorExpr
{
private:
  const FunctionSpaceExpr::expr_ptr trialSpace;
  const FunctionSpaceExpr::expr_ptr testSpace;
  const forms::BilinearFormIntegralSum sum;

  template<typename ForwardIterator>
  PropagationRules getPropagationRules(const ForwardIterator begin, const ForwardIterator end)
  {
    FormPropagationRulesGetter ruleGetter(*this);

    for(ForwardIterator i(begin); i!=end; ++i)
    {
      i->getTrialField()->accept(ruleGetter);
      i->getTestField()->accept(ruleGetter);
    }

    return ruleGetter.getRules();
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

  virtual PropagationRules getPropagationRules()
  {
    PropagationRules rules;
    rules += getPropagationRules(sum.begin_dx(), sum.end_dx());
    rules += getPropagationRules(sum.begin_ds(), sum.end_ds());
    rules += getPropagationRules(sum.begin_dS(), sum.end_dS());
    return rules;
  }
};

}

}

#endif
