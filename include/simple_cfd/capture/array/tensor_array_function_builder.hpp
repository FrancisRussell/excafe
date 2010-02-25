#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_BUILDER_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_BUILDER_HPP

#include <cstddef>
#include <map>
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
#include <simple_cfd/cell_manager.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class TensorArrayFunctionBuilderVisitor : public FieldVisitor
{
private:
  static const std::size_t dimension = D;
  typedef typename CellManager::ref<dimension>::general cell_ref_t;

  cell_ref_t cell;
  FreeTensorArray position;
  FreeTensorArray cellVertices;  
  std::map<Field, FreeTensorArray> discreteFieldCoefficients;
  std::stack<TensorFunction::ref> valueStack;

  TensorArrayFunctionPolynomial getVerticesAsTensorFunction() const
  {
    const ArrayIndex<fixed_tag> noArrayIndices(0);
    const TensorIndex<fixed_tag> noTensorIndices(0, dimension);

    TensorArrayFunctionPolynomial result(noArrayIndices, 0, dimension);
    result.appendVirtualTensorIndex();
    result.appendVirtualArrayIndex(cell->getCoordinateMapping().getSpaceDimension());

    const TensorFunction::polynomial_t poly = 
      ScalarReference(result.getIdentityArrayIndex(), result.getIdentityTensorIndex(), cellVertices);
    result(noArrayIndices, noTensorIndices) = poly;

    return result;
  }

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
  
  TensorArrayFunctionPolynomial buildCoordinateFunction() const
  {
    const FiniteElement<dimension>& coordinateMapping = cell->getCoordinateMapping();
    assert(coordinateMapping.getRank() == 0 && "Coordinate interpolation function must be rank 0");

    const TensorFunction::ref coordinateMappingBases(new TensorArrayFunctionPolynomial(coordinateMapping.getBasisFunctions(position)));
    const TensorFunction::ref coordinates(new TensorArrayFunctionPolynomial(getVerticesAsTensorFunction()));

    const std::size_t dimension = coordinateMappingBases->getTensorDimension();
    const std::size_t spaceDimension = coordinateMapping.getSpaceDimension();

    const ArrayIndex<fixed_tag> noArrayIndices(0);
    const TensorIndex<fixed_tag> noTensorIndices(0, dimension);

    TensorArrayFunctionProduct coordinateFunction(noArrayIndices, 0, dimension);
    const ArrayIndexID basisIndexID = coordinateFunction.newArrayIndex(spaceDimension);
    const TensorIndexID coordinateIndexID = coordinateFunction.newTensorIndex();

    const ArrayIndex<param_tag> basisIndex(1, &basisIndexID);
    const TensorIndex<param_tag> tensorIndex(1, &coordinateIndexID);
    coordinateFunction.addTerm(basisIndex, noTensorIndices, coordinateMappingBases);
    coordinateFunction.addTerm(basisIndex, tensorIndex, getVerticesAsTensorFunction());

    return coordinateFunction;
  }

  TensorFunction::ref buildInverseJacobian() const
  {
    TensorArrayFunctionPolynomial mapping(buildCoordinateFunction());
    //FIXME: implement me!
  }

public:
  TensorArrayFunctionBuilderVisitor(const cell_ref_t _cell, const FreeTensorArray& _position, 
    const FreeTensorArray& _cellVertices, 
    const std::map<Field, FreeTensorArray> _discreteFieldCoefficients) : 
    cell(_cell), position(_position), cellVertices(_cellVertices),
    discreteFieldCoefficients(_discreteFieldCoefficients)
  {
  }

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

    assert(arrayExtent == second->getArrayExtent());
    assert(rank == second->getTensorRank());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionSummation summation(arrayExtent, rank, dimension);
    const ArrayIndex<param_tag>& arrayIndex = summation.getIdentityArrayIndex();
    const TensorIndex<param_tag>& tensorIndex = summation.getIdentityTensorIndex();

    summation.addTerm(arrayIndex, tensorIndex, first);
    summation.addTerm(arrayIndex, tensorIndex, second);

    valueStack.push(TensorFunction::ref(new TensorArrayFunctionSummation(summation)));
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
    assert(arrayExtent == second->getArrayExtent());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionProduct product(arrayExtent, resultRank, dimension);
    const TensorIndexID summationIndex = product.newTensorIndex();

    TensorIndex<param_tag> firstTensorIndex = product.getIdentityTensorIndex().head(first->getTensorRank() - 1);
    firstTensorIndex.append(summationIndex);

    TensorIndex<param_tag> secondTensorIndex = product.getIdentityTensorIndex().tail(second->getTensorRank() - 1);
    secondTensorIndex.prepend(summationIndex);

    product.addTerm(product.getIdentityArrayIndex(), firstTensorIndex, first);
    product.addTerm(product.getIdentityArrayIndex(), secondTensorIndex, second);

    valueStack.push(TensorFunction::ref(new TensorArrayFunctionProduct(product)));
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
    assert(arrayExtent == second->getArrayExtent());
    assert(dimension == second->getTensorDimension());

    TensorArrayFunctionProduct product(arrayExtent, resultRank, dimension);

    const TensorIndex<param_tag> firstTensorIndex = product.getIdentityTensorIndex().head(first->getTensorRank());
    const TensorIndex<param_tag> secondTensorIndex = product.getIdentityTensorIndex().tail(second->getTensorRank());

    product.addTerm(product.getIdentityArrayIndex(), firstTensorIndex, first);
    product.addTerm(product.getIdentityArrayIndex(), secondTensorIndex, second);

    valueStack.push(TensorFunction::ref(new TensorArrayFunctionProduct(product)));
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
    assert(arrayExtent == second->getArrayExtent());
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

    valueStack.push(TensorFunction::ref(new TensorArrayFunctionProduct(product)));
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

  virtual void exit(FieldDivergence& divergence)
  {
    //FIXME: implement me!
  }

  // Terminals
  virtual void visit(FacetNormal& normal)
  {
    //FIXME: implement me!
  }

  virtual void visit(FieldBasis& basis)
  {
    //FIXME: implement me!
  }

  virtual void visit(FieldDiscreteReference& field)
  {
    //FIXME: implement me!
  }

  virtual void visit(FieldScalar& s)
  {
    //FIXME: implement me!
  }
};

}

}

#endif
