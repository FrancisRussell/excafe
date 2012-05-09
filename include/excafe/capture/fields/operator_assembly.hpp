#ifndef EXCAFE_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP
#define EXCAFE_CAPTURE_FIELDS_OPERATOR_ASSEMBLY_HPP

#include <set>
#include <memory>
#include "operator_expr.hpp"
#include "discrete_expr_visitor.hpp"
#include "function_space_expr.hpp"
#include <excafe/capture/forms/bilinear_form_integral_sum.hpp>
#include <excafe/capture/forms/field_discrete_reference_visitor.hpp>
#include <excafe/capture/indices/propagation_rules.hpp>

namespace excafe
{

namespace detail
{

class FormPropagationRulesGetter : public FieldDiscreteReferenceVisitor
{
private:
  OperatorAssembly& parent;
  PropagationRules rules;

public:
  FormPropagationRulesGetter(OperatorAssembly& _parent) : parent(_parent)
  {
  }
  
  PropagationRules getRules() const;
  void visit(FieldDiscreteReference& field);
  void visit(FieldScalar& s);
};

class FormDependencyGetter : public FieldDiscreteReferenceVisitor
{
private:
  std::set<DiscreteExpr*> dependencies;

public:
  std::set<DiscreteExpr*> getDependencies() const;
  void visit(FieldDiscreteReference& field);
  void visit(FieldScalar& s);
};


class OperatorAssembly : public OperatorExpr
{
private:
  const FunctionSpaceExpr::expr_ptr trialSpace;
  const FunctionSpaceExpr::expr_ptr testSpace;
  const forms::BilinearFormIntegralSum sum;

  template<typename ForwardIterator, typename Visitor>
  void applyFormVisitor(Visitor& v, const ForwardIterator begin, const ForwardIterator end) const
  {
    for(ForwardIterator i(begin); i!=end; ++i)
    {
      i->getTrialField()->accept(v);
      i->getTestField()->accept(v);
    }
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
    FormPropagationRulesGetter ruleGetter(*this);
    applyFormVisitor(ruleGetter, sum.begin_dx(), sum.end_dx());
    applyFormVisitor(ruleGetter, sum.begin_ds(), sum.end_ds());
    applyFormVisitor(ruleGetter, sum.begin_dS(), sum.end_dS());
    return ruleGetter.getRules();
  }

  virtual std::set<DiscreteExpr*> getDependencies() const 
  {
    FormDependencyGetter dependencyGetter;
    applyFormVisitor(dependencyGetter, sum.begin_dx(), sum.end_dx());
    applyFormVisitor(dependencyGetter, sum.begin_ds(), sum.end_ds());
    applyFormVisitor(dependencyGetter, sum.begin_dS(), sum.end_dS());
    return dependencyGetter.getDependencies();
  }
};

}

}

#endif
