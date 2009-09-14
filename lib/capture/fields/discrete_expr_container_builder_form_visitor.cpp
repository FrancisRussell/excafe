#include <simple_cfd/capture/fields/discrete_expr_container_builder_form_visitor.hpp>
#include <simple_cfd/capture/fields/discrete_expr_container_builder.hpp>
#include <simple_cfd/capture/forms/field_addition.hpp>
#include <simple_cfd/capture/forms/field_inner_product.hpp>
#include <simple_cfd/capture/forms/field_outer_product.hpp>
#include <simple_cfd/capture/forms/field_colon_product.hpp>
#include <simple_cfd/capture/forms/field_divergence.hpp>
#include <simple_cfd/capture/forms/field_gradient.hpp>

namespace cfd
{

namespace detail
{

DiscreteExprContainerBuilderFormVisitor::DiscreteExprContainerBuilderFormVisitor(DiscreteExprContainerBuilder& _parent) :
  parent(_parent)
{
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldAddition& addition)
{
  addition.getLeft().accept(*this);
  addition.getRight().accept(*this);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldInnerProduct& inner)
{
  inner.getLeft().accept(*this);
  inner.getRight().accept(*this);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldOuterProduct& outer)
{
  outer.getLeft().accept(*this);
  outer.getRight().accept(*this);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldColonProduct& colon)
{
  colon.getLeft().accept(*this);
  colon.getRight().accept(*this);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldGradient& gradient)
{
  gradient.getOperand().accept(*this);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldDivergence& divergence)
{
  divergence.getOperand().accept(*this);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FacetNormal& normal)
{
  // Nothing to do here
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldBasis& basis)
{
  // Nothing to do here
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldDiscreteReference& field) 
{
  field.getDiscreteField().getExpr()->accept(parent);
}

void DiscreteExprContainerBuilderFormVisitor::visit(FieldScalar& s)
{
  s.getValue().getExpr()->accept(parent);
}

}

}
