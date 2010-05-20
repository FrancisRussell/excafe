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
