#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_CONTAINER_BUILDER_FORM_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_CONTAINER_BUILDER_FORM_VISITOR_HPP

#include "fields_fwd.hpp"
#include <simple_cfd/capture/forms/forms_fwd.hpp>
#include <simple_cfd/capture/forms/field_visitor.hpp>

namespace cfd
{

namespace detail
{

class DiscreteExprContainerBuilderFormVisitor : public FieldVisitor
{
private:
  DiscreteExprContainerBuilder& parent;

public:
  DiscreteExprContainerBuilderFormVisitor(DiscreteExprContainerBuilder& _parent);
  virtual void visit(FieldAddition& addition);
  virtual void visit(FieldInnerProduct& inner);
  virtual void visit(FieldOuterProduct& outer);
  virtual void visit(FieldColonProduct& colon);
  virtual void visit(FieldGradient& gradient);
  virtual void visit(FieldDivergence& divergence);
  virtual void visit(FacetNormal& normal);
  virtual void visit(FieldBasis& basis);
  virtual void visit(FieldDiscreteReference& field);
  virtual void visit(FieldScalar& s);
};

}

}

#endif
