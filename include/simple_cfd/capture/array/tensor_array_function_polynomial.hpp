#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_POLYNOMIAL_HPP

#include <cstddef>
#include <functional>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <iterator>
#include <simple_cfd/exception.hpp>
#include <boost/foreach.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "array_traits.hpp"
#include "tensor_array_function.hpp"
#include "tensor_array_function_polynomial_visitor.hpp"
#include "tensor_array_function_helper.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionPolynomial : public TensorArrayFunction<TensorFunction::polynomial_t>
{
private:
  typedef TensorFunction::polynomial_t polynomial_t;
  typedef TensorArrayFunction<polynomial_t> parent_t;

  class ReferringIndicesCollector : public TensorArrayFunctionPolynomialVisitor
  {
  private:
    const FreeTensorArray freeTensorArray;
    std::set<ArrayIndexID> arrayIndices;
    std::set<TensorIndexID> tensorIndices;

  public:
    ReferringIndicesCollector(const FreeTensorArray& _freeTensorArray) : freeTensorArray(_freeTensorArray)
    {
    }
 
    void visit(const TensorArrayFunctionPolynomial& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      TensorFunction::polynomial_t& value)
    {
      const std::set<ScalarReference> references = value.getIndependentVariables();
      BOOST_FOREACH(const ScalarReference& reference, references)
      {
        if (!reference.isBound() && reference.getFreeTensorArray() == freeTensorArray && reference.isParameterised())
        {
          const ArrayIndex<param_tag> arrayIndex = reference.getArrayIndex();
          const std::set<ArrayIndexID> arrayParameters = arrayIndex.getReferencedParameters();
          arrayIndices.insert(arrayParameters.begin(), arrayParameters.end());
 
          const TensorIndex<param_tag> tensorIndex = reference.getTensorIndex();
          const std::set<TensorIndexID> tensorParameters = tensorIndex.getReferencedParameters();
          tensorIndices.insert(tensorIndices.begin(), tensorIndices.end());
        }
      }
    }
 
    std::set<ArrayIndexID> getArrayIndices() const
    {
      return arrayIndices;
    }
 
    std::set<TensorIndexID> getTensorIndices() const
    {
      return tensorIndices;
    }
  };

  class PolynomialCopierInternal : public TensorArrayFunctionPolynomialVisitor
  {
  private:
    const TensorArrayFunctionPolynomial& source;

  public:
    PolynomialCopierInternal(const TensorArrayFunctionPolynomial& _source) : source(_source)
    {
    }

    void visit(const TensorArrayFunctionPolynomial& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndex,
      const std::map<TensorIndexID, std::size_t>& tensorIndex,
      TensorFunction::polynomial_t& value)
    {
      value = source(arrayIndex, tensorIndex);
    }
  };

  class PolynomialCopierGeneral : public TensorArrayFunctionPolynomialVisitor
  {
  private:
    const TensorFunction& source;

  public:
    PolynomialCopierGeneral(const TensorFunction& _source) : source(_source)
    {
    }

    void visit(const TensorArrayFunctionPolynomial& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      TensorFunction::polynomial_t& value)
    {
      const ArrayIndex<fixed_tag> arrayIndex = 
        TensorArrayFunctionHelper::getIndex(arrayIndexMap, parent.getIdentityArrayIndex());
      const TensorIndex<fixed_tag> tensorIndex = 
        TensorArrayFunctionHelper::getIndex(tensorIndexMap, parent.getIdentityTensorIndex());

      value = source.getPolynomial(arrayIndex, tensorIndex);
    }
  };


  class IndexSpecialiser : public TensorArrayFunctionPolynomialVisitor
  {
  public:
    void visit(const TensorArrayFunctionPolynomial& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      TensorFunction::polynomial_t& value)
    {
      const std::set<ScalarReference> references = value.getIndependentVariables();
      BOOST_FOREACH(const ScalarReference& reference, references)
      {
        if (reference.isParameterised())
        {
          const ArrayIndex<param_tag> arrayIndex = reference.getArrayIndex();
          const ArrayIndex<param_tag> newArrayIndex = arrayIndex.substituteLiterals(arrayIndexMap);
 
          const TensorIndex<param_tag> tensorIndex = reference.getTensorIndex();
          const TensorIndex<param_tag> newTensorIndex = tensorIndex.substituteLiterals(tensorIndexMap);
          
          const ScalarReference newReference(arrayIndex, tensorIndex, reference);
          value.replaceIndependentVariable(reference, newReference);
        }
      }
    }
  };

  class Differentiator : public TensorArrayFunctionPolynomialVisitor
  {
  private:
    const ScalarReference variable;

  public:
    Differentiator(const ScalarReference& _variable) : variable(_variable)
    {
      assert(!variable.isBound() && !variable.isParameterised());
    }

    void visit(const TensorArrayFunctionPolynomial& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      TensorFunction::polynomial_t& value)
    {
      value = value.derivative(variable);
    }
  };

  TensorArrayFunctionPolynomial(const TensorArrayFunctionPolynomial& original,
    const std::set<ArrayIndexID>& expandArrayIndices, const std::set<TensorIndexID>& expandTensorIndices) :
    parent_t(original.getArrayExtent(), original.getTensorRank(), original.getTensorDimension(),
    original.arrayIndexParameters, original.tensorIndexParameters)
  {
    // Calculate virtual parameters
    std::set_difference(original.arrayVirtualParameters.begin(), original.arrayVirtualParameters.end(),
      expandArrayIndices.begin(), expandArrayIndices.end(),
      std::inserter(arrayVirtualParameters, arrayVirtualParameters.begin()));

    std::set_difference(original.tensorVirtualParameters.begin(), original.tensorVirtualParameters.end(),
      expandTensorIndices.begin(), expandTensorIndices.end(),
      std::inserter(tensorVirtualParameters, tensorVirtualParameters.begin()));

    // Resize values vector to correct extent
    values.resize(extentReal());

    // Use visitor to perform polynomial assignments
    PolynomialCopierInternal copier(original);
    accept(copier);

    // Specialise known indices
    specialiseKnownIndices();
  }

  void specialiseKnownIndices() 
  {
    IndexSpecialiser specialiser;
    accept(specialiser);
  }

public:
  TensorArrayFunctionPolynomial(const ArrayIndex<fixed_tag>& _arrayExtents, const std::size_t _rank, 
    const std::size_t _dimension) : parent_t(_arrayExtents, _rank, _dimension)
  {
  }

  TensorArrayFunctionPolynomial(const TensorFunction& original) : 
    parent_t(original.getArrayExtent(), original.getTensorRank(), original.getTensorDimension())
  {
    // Use visitor to perform polynomial assignments
    PolynomialCopierGeneral copier(original);
    accept(copier);

    // Specialise known indices
    specialiseKnownIndices();
  }

  void accept(TensorArrayFunctionPolynomialVisitor& visitor)
  {
    parent_t::IndexIncrementer incrementer(*this);
    std::map<ArrayIndexID, std::size_t> arrayIndex; 
    std::map<TensorIndexID, std::size_t> tensorIndex;

    incrementer.zero(arrayIndex);
    incrementer.zero(tensorIndex);

    do
    {
      visitor.visit(*this, arrayIndex, tensorIndex, (*this)(arrayIndex, tensorIndex));
    }
    while(!incrementer.increment(arrayIndex, tensorIndex));
  }

  virtual polynomial_t getPolynomial(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::map<ArrayIndexID, std::size_t> arrayIndexMap =
      TensorArrayFunctionHelper::indexToMap(arrayIndexParameters, arrayIndex);

    const std::map<TensorIndexID, std::size_t> tensorIndexMap =
      TensorArrayFunctionHelper::indexToMap(tensorIndexParameters, tensorIndex);

    polynomial_t poly = (*this)(arrayIndex, tensorIndex);
    IndexSpecialiser specialiser;
    specialiser.visit(*this, arrayIndexMap, tensorIndexMap, poly);

    //TODO: assert that this polynomial isn't parameterised
    return poly;
  }

  virtual TensorFunction::ref differentiate(const ScalarReference& reference)
  {
    if (reference.isBound()) 
        CFD_EXCEPTION("Cannot differentiate with respect to bound reference.");

    if (reference.isParameterised()) 
      CFD_EXCEPTION("Cannot differentiate with respect to parameterised reference.");

    ReferringIndicesCollector indexCollector(reference.getFreeTensorArray());
    accept(indexCollector);

    TensorArrayFunctionPolynomial result(*this, indexCollector.getArrayIndices(),
      indexCollector.getTensorIndices());

    Differentiator differentiator(reference);
    result.accept(differentiator);

    return TensorFunction::ref(new TensorArrayFunctionPolynomial(result));
  }
};

}

}

#endif
