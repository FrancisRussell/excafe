#ifndef SIMPLE_CFD_FORMS_FIELD_VISITOR_HPP
#define SIMPLE_CFD_FORMS_FIELD_VISITOR_HPP

#include "forms_fwd.hpp"

namespace cfd
{

namespace forms
{

class FieldVisitor
{
public:
  // Non Terminals
  virtual void enter(Addition& addition) = 0;
  virtual void exit(Addition& addition) = 0;

  virtual void enter(InnerProduct& inner) = 0;
  virtual void exit(InnerProduct& inner) = 0;

  virtual void enter(OuterProduct& outer) = 0;
  virtual void exit(OuterProduct& outer) = 0;

  virtual void enter(ColonProduct& colon) = 0;
  virtual void exit(ColonProduct& colon) = 0;

  virtual void enter(Gradient& gradient) = 0;
  virtual void exit(Gradient& gradient) = 0;

  virtual void enter(Divergence& divergence) = 0;
  virtual void exit(Divergence& divergence) = 0;

  // Terminals
  virtual void visit(FacetNormal& normal) = 0;
  virtual void visit(BasisField& basis) = 0;
  virtual void visit(DiscreteFieldReference& field) = 0;
  virtual void visit(TensorLiteral& literal) = 0;

  virtual ~FieldVisitor() {}
};

}

}

#endif
