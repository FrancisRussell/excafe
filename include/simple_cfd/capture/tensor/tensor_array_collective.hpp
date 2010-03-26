#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_COLLECTIVE_HPP

#include <vector>
#include <cassert>
#include "tensor_fwd.hpp"
#include "tensor_array.hpp"
#include "tensor_array_ref.hpp"
#include "index.hpp"
#include "indexed_tensor.hpp"
#include "index_generator.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayCollective : public TensorArray
{
private:
  IndexGenerator generator;
  TensorIndex visibleIndices;
  TensorIndex hiddenIndices;
  std::vector<IndexedTensor> operands;

  std::size_t getDimension() const
  {
    return getTensorSize().getDimension();
  }

  void generateVisibleTensorIndices(IndexGenerator& g)
  {
    const TensorSize tensorSize = visibleIndices.getSize();
    for(std::size_t i=0; i<tensorSize.numIndices(); ++i)
    {
      visibleIndices[i] = g.newTensorIndexVariable(tensorSize.getLimit(i));
    }
  }

protected:
  std::size_t numOperands() const
  {
    return operands.size();
  }

  IndexedTensor& operand(const std::size_t n)
  {
    return operands[n];
  }

  const IndexedTensor& operand(const std::size_t n) const
  {
    return operands[n];
  }

public:
  typedef std::vector<IndexedTensor>::iterator iterator;
  typedef std::vector<IndexedTensor>::const_iterator const_iterator;

  TensorArrayCollective(const TensorSize& tensorSize) :
    visibleIndices(tensorSize), hiddenIndices(TensorSize(0, getDimension()))
  {
    generateVisibleTensorIndices(generator);
    assert(visibleIndices.allVariable());
  }

  TensorArrayCollective(const TensorArrayCollective& t) : generator(t.generator),
    visibleIndices(t.visibleIndices), hiddenIndices(t.hiddenIndices), 
    operands(t.operands)
  {
  }

  iterator begin()
  {
    return operands.begin();
  }

  const_iterator begin() const
  {
    return operands.begin();
  }

  iterator end()
  {
    return operands.end();
  }

  const_iterator end() const
  {
    return operands.end();
  }

  TensorSize getTensorSize() const
  {
    return visibleIndices.getSize();
  }

  TensorIndexVariable newHiddenIndex()
  {
    const TensorIndexVariable newIndex = generator.newTensorIndexVariable(getDimension());
    hiddenIndices = hiddenIndices.append(newIndex);
    return newIndex;
  }

  void createNewHiddenIndices(const std::size_t n) 
  {
    for(std::size_t i=0; i<n; ++i)
      newHiddenIndex();
  }

  TensorIndex getVisibleIndices() const
  {
    return visibleIndices;
  }

  TensorIndex getHiddenIndices() const
  {
    return hiddenIndices;
  }

  void addOperand(const TensorIndex& index, const TensorArrayRef& tensor)
  {
    //TODO: check index only uses indices we've defined somewhere
    operands.push_back(IndexedTensor(tensor, index));
  }
};

}

}

#endif
