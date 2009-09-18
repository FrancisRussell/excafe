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

  virtual void enter(FieldAddition& addition);
  virtual void exit(FieldAddition& addition);

  virtual void enter(FieldInnerProduct& inner);
  virtual void exit(FieldInnerProduct& inner);

  virtual void enter(FieldOuterProduct& outer);
  virtual void exit(FieldOuterProduct& outer);

  virtual void enter(FieldColonProduct& colon);
  virtual void exit(FieldColonProduct& colon);

  virtual void enter(FieldGradient& gradient);
  virtual void exit(FieldGradient& gradient);

  virtual void enter(FieldDivergence& divergence);
  virtual void exit(FieldDivergence& divergence);

  virtual void visit(FacetNormal& normal);
  virtual void visit(FieldBasis& basis);
  virtual void visit(FieldDiscreteReference& field);
  virtual void visit(FieldScalar& s);
};

}

}

#endif
