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
#include <boost/operators.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "array_traits.hpp"
#include "tensor_array_function.hpp"
#include "tensor_array_function_visitor.hpp"
#include "tensor_array_function_helper.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayFunctionPolynomial : public TensorArrayFunction<TensorFunction::polynomial_t>, 
                                      boost::multipliable<TensorArrayFunctionPolynomial,
                                      boost::multipliable<TensorArrayFunctionPolynomial, double,
                                      boost::addable<TensorArrayFunctionPolynomial
                                      > > >
{
private:
  typedef element_t polynomial_t;
  typedef TensorArrayFunction<polynomial_t> parent_t;
  typedef TensorArrayFunctionVisitor<polynomial_t> visitor_t;

  class ReferringIndicesCollector : public visitor_t
  {
  private:
    const FreeTensorArray freeTensorArray;
    std::set<ArrayIndexID> arrayIndices;
    std::set<TensorIndexID> tensorIndices;

  public:
    ReferringIndicesCollector(const FreeTensorArray& _freeTensorArray) : freeTensorArray(_freeTensorArray)
    {
    }
 
    void visit(const TensorArrayFunction<element_t>& parent,
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

  class PolynomialCopierInternal : public visitor_t
  {
  private:
    const TensorArrayFunctionPolynomial& source;

  public:
    PolynomialCopierInternal(const TensorArrayFunctionPolynomial& _source) : source(_source)
    {
    }

    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndex,
      const std::map<TensorIndexID, std::size_t>& tensorIndex,
      polynomial_t& value)
    {
      value = source(arrayIndex, tensorIndex);
    }
  };

  class PolynomialCopierGeneral : public visitor_t
  {
  private:
    const TensorFunction& source;

  public:
    PolynomialCopierGeneral(const TensorFunction& _source) : source(_source)
    {
    }

    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      polynomial_t& value)
    {
      const ArrayIndex<fixed_tag> arrayIndex = 
        TensorArrayFunctionHelper::getIndex(arrayIndexMap, parent.getIdentityArrayIndex());
      const TensorIndex<fixed_tag> tensorIndex = 
        TensorArrayFunctionHelper::getIndex(tensorIndexMap, parent.getIdentityTensorIndex());

      value = source.getPolynomial(arrayIndex, tensorIndex);
    }
  };

  class SingleTensorEntryCopier: public visitor_t
  {
  private:
    const TensorFunction& source;
    const TensorIndex<fixed_tag> sourceIndex;

  public:
    SingleTensorEntryCopier(const TensorFunction& _source, const TensorIndex<fixed_tag>& _sourceIndex) : 
    source(_source), sourceIndex(_sourceIndex)
    {
    }

    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      polynomial_t& value)
    {
      const ArrayIndex<fixed_tag> arrayIndex = 
        TensorArrayFunctionHelper::getIndex(arrayIndexMap, parent.getIdentityArrayIndex());

      value = source.getPolynomial(arrayIndex, sourceIndex);
    }
  };

  class IndexSpecialiser : public visitor_t
  {
  public:
    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      polynomial_t& value)
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

  class Differentiator : public visitor_t
  {
  private:
    const ScalarReference variable;

  public:
    Differentiator(const ScalarReference& _variable) : variable(_variable)
    {
      assert(!variable.isBound() && !variable.isParameterised());
    }

    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      polynomial_t& value)
    {
      value = value.derivative(variable);
    }
  };

  class SubMatrixFiller : public visitor_t
  {
  private:
    const TensorFunction& source;
    const std::size_t ignoredRow;
    const std::size_t ignoredCol;

  public:
    SubMatrixFiller(const TensorFunction& _source, const std::size_t row, const std::size_t col) : 
      source(_source), ignoredRow(row), ignoredCol(col)
    {
      assert(source.getTensorRank() == 2);
    }

    void visit(const TensorArrayFunction<element_t>& parent,
      const std::map<ArrayIndexID, std::size_t>& arrayIndexMap,
      const std::map<TensorIndexID, std::size_t>& tensorIndexMap,
      polynomial_t& value)
    {
      const ArrayIndex<fixed_tag> arrayIndex = 
        TensorArrayFunctionHelper::getIndex(arrayIndexMap, parent.getIdentityArrayIndex());
      const TensorIndex<fixed_tag> tensorIndex = 
        TensorArrayFunctionHelper::getIndex(tensorIndexMap, parent.getIdentityTensorIndex());

      TensorIndex<fixed_tag> modifiedTensorIndex(tensorIndex);
      if (modifiedTensorIndex[0] >= ignoredRow) ++modifiedTensorIndex[0];
      if (modifiedTensorIndex[1] >= ignoredCol) ++modifiedTensorIndex[1];

      value = source.getPolynomial(arrayIndex, modifiedTensorIndex);
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

  virtual polynomial_t getPolynomial(const ArrayIndex<fixed_tag>& arrayIndex, const TensorIndex<fixed_tag>& tensorIndex) const
  {
    const std::map<ArrayIndexID, std::size_t> arrayIndexMap =
      TensorArrayFunctionHelper::indexToMap(arrayIndexParameters, arrayIndex);

    const std::map<TensorIndexID, std::size_t> tensorIndexMap =
      TensorArrayFunctionHelper::indexToMap(tensorIndexParameters, tensorIndex);

    polynomial_t poly = (*this)(arrayIndexMap, tensorIndexMap);
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
  
  TensorArrayFunctionPolynomial& operator*=(const double s)
  {
    BOOST_FOREACH(polynomial_t& p, values)
    {
      p *= s;
    }
    return *this;
  }

  TensorArrayFunctionPolynomial& operator+=(const TensorFunction& f)
  {
    assert(f.getArrayExtent() == getArrayExtent());
    assert(f.getTensorRank() == getTensorRank());
    assert(f.getTensorDimension() == getTensorDimension());

    // Double check that this is really doing an expansion
    TensorArrayFunctionPolynomial expandedThis(*this);
    const TensorArrayFunctionPolynomial expandedF(f);

    std::transform(expandedThis.values.begin(), expandedThis.values.end(), expandedF.values.begin(),
      expandedThis.values.begin(), std::plus<polynomial_t>());

    std::swap(*this, expandedThis);
    return *this;
  }

  TensorArrayFunctionPolynomial& operator*=(const TensorFunction& f)
  {
    assert(f.getArrayExtent() == getArrayExtent());
    assert(f.getTensorRank() == getTensorRank());
    assert(f.getTensorDimension() == getTensorDimension());

    // Double check that this is really doing an expansion
    TensorArrayFunctionPolynomial expandedThis(*this);
    const TensorArrayFunctionPolynomial expandedF(f);

    std::transform(expandedThis.values.begin(), expandedThis.values.end(), expandedF.values.begin(),
      expandedThis.values.begin(), std::multiplies<polynomial_t>());

    std::swap(*this, expandedThis);
    return *this;
  }

  TensorArrayFunctionPolynomial getSubMatrices(const std::size_t row, const std::size_t col) const
  {
    assert(getTensorDimension() >= 2);
    assert(getTensorRank() == 2);
    assert(row < dimension);
    assert(col < dimension);

    const ArrayIndex<fixed_tag> arrayExtents = getArrayExtent();

    TensorArrayFunctionPolynomial result(arrayExtents, getTensorRank(), getTensorDimension()-1);
    SubMatrixFiller subMatrixFiller(*this, row, col);
    result.accept(subMatrixFiller);

    return result;
  }

  TensorArrayFunctionPolynomial getTensorEntry(const TensorIndex<fixed_tag>& index) const
  {
    assert(index.getRank() == getTensorRank());
    assert(index.getDimension() == getTensorDimension());

    TensorArrayFunctionPolynomial result(getArrayExtent(), 0, dimension);
    SingleTensorEntryCopier copier(*this, index);
    result.accept(copier);
    return result;
  }

  TensorArrayFunctionPolynomial getDeterminants() const
  {
    assert(getTensorRank() == 2);

    const std::size_t dimension = getTensorDimension();
    const ArrayIndex<fixed_tag> arrayExtents = getArrayExtent();

    if (dimension == 1)
    {
      return *this;
    }
    else
    {
      TensorArrayFunctionPolynomial result(arrayExtents, 0, dimension);

      for(std::size_t i=0; i<dimension; ++i)
      {
        for(std::size_t j=0; j<dimension; ++i)
        {
          TensorIndex<fixed_tag> tensorIndex(2, dimension);
          tensorIndex[0] = i;
          tensorIndex[1] = j;

          const double sign = (i+j)%2 == 0 ? 1.0 : -1.0;
          result += sign * getTensorEntry(tensorIndex) * getMinors(i, j);
        }
      }

      return result;
    }
  }

  TensorArrayFunctionPolynomial getAdjugateMatrix() const
  {
    assert(getTensorRank() == 2);

    const std::size_t dimension = getTensorDimension();
    const ArrayIndex<fixed_tag> arrayExtents = getArrayExtent();
    const std::vector<std::size_t> arrayExtentsVector(arrayExtents.begin(), arrayExtents.end());

    const TensorIndex<fixed_tag> noTensorIndex(0, dimension);
    const std::vector<TensorIndexID> noTensorIndices;

    TensorArrayFunctionPolynomial result(arrayExtents, 2, dimension);

    for(std::size_t i=0; i<dimension; ++i)
    {
      for(std::size_t j=0; j<dimension; ++i)
      {
        TensorIndex<fixed_tag> currentTensorIndex(2, dimension);
        // Transpose of adjugate
        currentTensorIndex[0]= j;
        currentTensorIndex[1]= i;

        const double sign = (i+j)%2 == 0 ? 1.0 : -1.0;
        const TensorArrayFunctionPolynomial minors = getMinors(i, j);

        IndexIncrementer incrementer(arrayIndexParameters, noTensorIndices, arrayExtentsVector, dimension);

        std::map<ArrayIndexID, std::size_t> arrayIndexMap;
        std::map<TensorIndexID, std::size_t> tensorIndexMap;

        incrementer.zero(arrayIndexMap);
        incrementer.zero(tensorIndexMap);

        do
        {
          const ArrayIndex<fixed_tag> arrayIndex = TensorArrayFunctionHelper::getIndex(arrayIndexMap,
            getIdentityArrayIndex());
          result(arrayIndex, currentTensorIndex) += sign * minors(arrayIndex, noTensorIndex);
        }
        while(!incrementer.increment(arrayIndexMap, tensorIndexMap));
      }
    }

    return result;
  }
 

  TensorArrayFunctionPolynomial getMinors(const std::size_t row, const std::size_t col) const
  {
     const TensorArrayFunctionPolynomial subMatrices = getSubMatrices(row, col);
     const TensorArrayFunctionPolynomial determinants = subMatrices.getDeterminants();
     return determinants;
  }
};

}

}

#endif
