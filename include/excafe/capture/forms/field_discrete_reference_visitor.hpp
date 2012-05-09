#ifndef EXCAFE_FORMS_FIELD_DISCRETE_REFERENCE_VISITOR_HPP
#define EXCAFE_FORMS_FIELD_DISCRETE_REFERENCE_VISITOR_HPP

#include "forms_fwd.hpp"
#include "field_visitor.hpp"

namespace excafe
{

namespace detail
{

class FieldDiscreteReferenceVisitor : public FieldVisitor
{
public:
  // Non Terminals
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
};

}

}



#endif
