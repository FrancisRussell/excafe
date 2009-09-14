#ifndef SIMPLE_CFD_FORMS_FIELD_VISITOR_HPP
#define SIMPLE_CFD_FORMS_FIELD_VISITOR_HPP

#include "forms_fwd.hpp"

namespace cfd
{

namespace detail
{

class FieldVisitor
{
public:
  // Non Terminals
  virtual void visit(FieldAddition& addition) = 0;
  virtual void visit(FieldInnerProduct& inner) = 0;
  virtual void visit(FieldOuterProduct& outer) = 0;
  virtual void visit(FieldColonProduct& colon) = 0;
  virtual void visit(FieldGradient& gradient) = 0;
  virtual void visit(FieldDivergence& divergence) = 0;

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
