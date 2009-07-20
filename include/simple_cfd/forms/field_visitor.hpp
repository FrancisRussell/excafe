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

  virtual void enter(InnerProduct& addition) = 0;
  virtual void exit(InnerProduct& addition) = 0;

  virtual void enter(OuterProduct& addition) = 0;
  virtual void exit(OuterProduct& addition) = 0;

  virtual void enter(ColonProduct& addition) = 0;
  virtual void exit(ColonProduct& addition) = 0;

  virtual void enter(Gradient& gradient) = 0;
  virtual void exit(Gradient& gradient) = 0;

  virtual void enter(Divergence& gradient) = 0;
  virtual void exit(Divergence& gradient) = 0;

  // Terminals
  virtual void visit(BasisField& basis) = 0;
  virtual void visit(DiscreteField& basis) = 0;
  virtual void visit(TensorLiteral& basis) = 0;

  virtual ~FieldVisitor();
};

}

}

#endif
