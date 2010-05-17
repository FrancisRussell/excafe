#ifndef SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_CONTAINER_BUILDER_FORM_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_DISCRETE_EXPR_CONTAINER_BUILDER_FORM_VISITOR_HPP

#include "fields_fwd.hpp"
#include <simple_cfd/capture/forms/forms_fwd.hpp>
#include <simple_cfd/capture/forms/field_discrete_reference_visitor.hpp>

namespace cfd
{

namespace detail
{

class DiscreteExprContainerBuilderFormVisitor : public FieldDiscreteReferenceVisitor
{
private:
  DiscreteExprContainerBuilder& parent;

public:
  DiscreteExprContainerBuilderFormVisitor(DiscreteExprContainerBuilder& _parent);
  virtual void visit(FieldDiscreteReference& field);
  virtual void visit(FieldScalar& s);
};

}

}

#endif
