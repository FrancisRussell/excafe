#include <memory>
#include <set>
#include <simple_cfd/capture/fields/operator_assembly.hpp>
#include <simple_cfd/capture/indices/index_propagation_all.hpp>

namespace cfd
{

namespace detail
{

PropagationRules FormPropagationRulesGetter::getRules() const
{
  return rules;
}

void FormPropagationRulesGetter::visit(FieldDiscreteReference& field)
{
  rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*field.getDiscreteField().getExpr(), parent)));
}

void FormPropagationRulesGetter::visit(FieldScalar& s)
{
  rules.insert(std::auto_ptr<PropagationRule>(new IndexPropagationAll(*s.getValue().getExpr(), parent)));
}


std::set<DiscreteExpr*> FormDependencyGetter::getDependencies() const
{
  return dependencies;
}

void FormDependencyGetter::visit(FieldDiscreteReference& field)
{
  dependencies.insert(&(*field.getDiscreteField().getExpr()));
}

void FormDependencyGetter::visit(FieldScalar& s)
{
  dependencies.insert(&(*s.getValue().getExpr()));
}

}

}
