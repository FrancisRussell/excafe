#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_BUILDER
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_BUILDER

#include <cstddef>
#include <map>
#include <memory>
#include "free_tensor_array.hpp"
#include <simple_cfd/capture/fields/field.hpp>
#include <simple_cfd/capture/forms/field_visitor.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class TensorArrayFunctionBuilderVisitor : public FieldVisitor
{
private:
  static const std::size_t dimension = D;

  const std::auto_ptr< GeneralCell<dimension> > cell;
  FreeTensorArray position;
  FreeTensorArray cellVertices;  
  std::map<Field, FreeTensorArray> discreteFieldCoefficients;

public:
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
};

}

}

#endif
