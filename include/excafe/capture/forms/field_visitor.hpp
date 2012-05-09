#ifndef EXCAFE_FORMS_FIELD_VISITOR_HPP
#define EXCAFE_FORMS_FIELD_VISITOR_HPP

#include "forms_fwd.hpp"

namespace excafe
{

namespace detail
{

class FieldVisitor
{
public:
  // Non Terminals
  virtual void enter(FieldAddition& addition) = 0;
  virtual void exit(FieldAddition& addition) = 0;

  virtual void enter(FieldInnerProduct& inner) = 0;
  virtual void exit(FieldInnerProduct& inner) = 0;

  virtual void enter(FieldOuterProduct& outer) = 0;
  virtual void exit(FieldOuterProduct& outer) = 0;

  virtual void enter(FieldColonProduct& colon) = 0;
  virtual void exit(FieldColonProduct& colon) = 0;

  virtual void enter(FieldGradient& gradient) = 0;
  virtual void exit(FieldGradient& gradient) = 0;

  virtual void enter(FieldDivergence& divergence) = 0;
  virtual void exit(FieldDivergence& divergence) = 0;

  // Terminals
  virtual void visit(FacetNormal& normal) = 0;
  virtual void visit(FieldBasis& basis) = 0;
  virtual void visit(FieldDiscreteReference& field) = 0;
  virtual void visit(FieldScalar& s) = 0;

  virtual ~FieldVisitor() {}
};

}

}

#endif
