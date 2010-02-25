#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_BUILDER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_BUILDER_HPP

#include <cstddef>
#include <map>
#include <memory>
#include <stack>
#include <cassert>
#include "free_tensor_array.hpp"
#include "tensor_array_function_references.hpp"
#include "tensor_array_function_call.hpp"
#include "scalar_reference.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
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
  std::stack<TensorFunction::ref> valueStack;

  TensorFunction::ref buildGradient(const TensorFunction::ref function) const
  {
    const ArrayIndex<fixed_tag> arrayExtents = function->getArrayExtent();
    const std::size_t rank = function->getTensorRank();
    const std::size_t dimension = function->getTensorDimension();

    const ArrayIndex<fixed_tag> noArrayIndices(0);

    TensorArrayFunctionReferences gradientTensor(noArrayIndices, 1, dimension);
    const std::vector<ArrayIndexID> virtualArrayIndices = gradientTensor.appendVirtualArrayIndices(arrayExtents);
    const std::vector<TensorIndexID> virtualTensorIndices = gradientTensor.appendVirtualTensorIndices(rank);

    const ArrayIndex<fixed_tag> virtualArrayIndex(virtualArrayIndices.size(), &virtualArrayIndices[0]);
    const TensorIndex<fixed_tag> virtualTensorIndex(rank, dimension, &virtualTensorIndices[0]);

    for(std::size_t d=0; d<dimension; ++d)
    {
      const TensorIndex<fixed_tag> realTensorIndex(1, dimension, &d);
      const ScalarReference coordinate(noArrayIndices, realTensorIndex, position);
      gradientTensor(noArrayIndices, realTensorIndex) = 
        TensorArrayFunctionCall(virtualArrayIndex, virtualTensorIndex, function->differentiate(coordinate));
    }

    return TensorFunction::ref(new TensorArrayFunctionReferences(gradientTensor));
  }
  
  TensorArrayFunctionPolynomial buildJacobian() const
  {
  }

public:
  virtual void enter(FieldAddition& addition)
  {
  }

  virtual void exit(FieldAddition& addition)
  {
    const TensorFunction::ref second = valueStack.top();
    valueStack.pop();

    const TensorFunction::ref first = valueStack.top();
    valueStack.pop();

    const ArrayIndex<fixed_tag> arrayExtent = first->getArrayExtent();
    const std::size_t rank = first->getTensorRank();
    const std::size_t dimension = first->getTensorRank();

    assert(arrayExtent = second->getArrayExtent());
    assert(rank == second->getTensorRank());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionSummation summation(arrayExtent, rank, dimension);
    const ArrayIndex<param_tag>& arrayIndex = summation.getIdentityArrayIndex();
    const TensorIndex<param_tag>& tensorIndex = summation.getIdentityTensorIndex();

    summation.addTerm(arrayIndex, tensorIndex, first);
    summation.addTerm(arrayIndex, tensorIndex, second);

    return TensorFunction::ref(new TensorArrayFunctionSummation(summation));
  }

  virtual void enter(FieldInnerProduct& inner)
  {
  }

  virtual void exit(FieldInnerProduct& inner)
  {
    const TensorFunction::ref second = valueStack.top();
    valueStack.pop();

    const TensorFunction::ref first = valueStack.top();
    valueStack.pop();

    const ArrayIndex<fixed_tag> arrayExtent = first->getArrayExtent();
    const std::size_t dimension = first->getTensorDimension();
    const std::size_t resultRank = first->getTensorRank() + second->getTensorRank() - 2;

    assert(arrayExtent.numIndices() == 1);
    assert(arrayExtent = second->getArrayExtent());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionProduct product(arrayExtent, resultRank, dimension);
    const TensorIndexID summationIndex = product.newTensorIndex();

    TensorIndex<param_tag> firstTensorIndex = product.getIdentityTensorIndex().head(first->getTensorRank() - 1);
    firstTensorIndex.append(summationIndex);

    TensorIndex<param_tag> secondTensorIndex = product.getIdentityTensorIndex().tail(second->getTensorRank() - 1);
    secondTensorIndex.prepend(summationIndex);

    product.addTerm(product.getIdentityArrayIndex(), firstTensorIndex, first);
    product.addTerm(product.getIdentityArrayIndex(), secondTensorIndex, second);

    return TensorFunction::ref(new TensorArrayFunctionProduct(product));
  }

  virtual void enter(FieldOuterProduct& outer)
  {
  }

  virtual void exit(FieldOuterProduct& outer)
  {
    const TensorFunction::ref second = valueStack.top();
    valueStack.pop();

    const TensorFunction::ref first = valueStack.top();
    valueStack.pop();

    const ArrayIndex<fixed_tag> arrayExtent = first->getArrayExtent();
    const std::size_t dimension = first->getTensorDimension();
    const std::size_t resultRank = first->getTensorRank() + second->getTensorRank();

    assert(arrayExtent.numIndices() == 1);
    assert(arrayExtent = second->getArrayExtent());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionProduct product(arrayExtent, resultRank, dimension);

    const TensorIndex<param_tag> firstTensorIndex = product.getIdentityTensorIndex().head(first->getTensorRank());
    const TensorIndex<param_tag> secondTensorIndex = product.getIdentityTensorIndex().tail(second->getTensorRank());

    product.addTerm(product.getIdentityArrayIndex(), firstTensorIndex, first);
    product.addTerm(product.getIdentityArrayIndex(), secondTensorIndex, second);

    return TensorFunction::ref(new TensorArrayFunctionProduct(product));
  }

  virtual void enter(FieldColonProduct& colon)
  {
  }

  virtual void exit(FieldColonProduct& colon)
  {
    const TensorFunction::ref second = valueStack.top();
    valueStack.pop();

    const TensorFunction::ref first = valueStack.top();
    valueStack.pop();

    const ArrayIndex<fixed_tag> arrayExtent = first->getArrayExtent();
    const std::size_t dimension = first->getTensorDimension();
    const std::size_t resultRank = 0;

    assert(arrayExtent.numIndices() == 1);
    assert(arrayExtent = second->getArrayExtent());
    assert(first->getTensorRank() == second->getTensorRank());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionProduct product(arrayExtent, resultRank, dimension);

    const std::vector<TensorIndexID> firstTensorParams = product.newTensorIndices(first->getTensorRank());
    const std::vector<TensorIndexID> secondTensorParams = product.newTensorIndices(second->getTensorRank());

    const TensorIndex<param_tag> firstTensorIndex = 
      TensorIndex<param_tag>(firstTensorParams.size(), dimension, &firstTensorParams[0]);

    const TensorIndex<param_tag> secondTensorIndex = 
      TensorIndex<param_tag>(secondTensorParams.size(), dimension, &secondTensorParams[0]);

    product.addTerm(product.getIdentityArrayIndex(), firstTensorIndex, first);
    product.addTerm(product.getIdentityArrayIndex(), secondTensorIndex, second);

    return TensorFunction::ref(new TensorArrayFunctionProduct(product));
  }

  virtual void enter(FieldGradient& gradient)
  {
  }

  virtual void exit(FieldGradient& gradient)
  {
    //FIXME: implement me!

    const TensorFunction::ref value = valueStack.top();
    valueStack.pop();
    
    const TensorFunction::ref localGradient = buildGradient(value);
  }

  virtual void enter(FieldDivergence& divergence)
  {
  }

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
