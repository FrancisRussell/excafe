#ifndef SIMPLE_CFD_NUMERIC_TENSOR_HPP
#define SIMPLE_CFD_NUMERIC_TENSOR_HPP

#include <simple_cfd_fwd.hpp>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <functional>
#include <numeric>
#include <algorithm>
#include <functional>
#include <ostream>
#include "index.hpp"
#include "tensor_size.hpp"

namespace cfd
{

namespace detail
{

template<unsigned  X, unsigned Y>
struct Power
{
  static const unsigned int value = X * Power<X, Y-1>::value;
};

template<unsigned X>
struct Power<X, 0>
{
  static const unsigned int value = 1;
};

}

template<std::size_t D, typename T>
class Tensor
{
public:
  static const std::size_t  dimension = D;
  typedef std::size_t size_type;
  typedef T value_type;
  typedef row_major_tag layout_tag;
  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::iterator const_iterator;

  TensorSize size;

private:
  std::vector<value_type>  elements;

  static std::size_t pow(const std::size_t base, const std::size_t exponent)
  {
    std::size_t result = 1;

    for(std::size_t i=0; i<exponent; ++i)
      result *= base;

    return result;
  }

public:
  Tensor() : size(0, dimension),  elements(getExtent())
  {
    std::fill(elements.begin(), elements.end(), value_type());
  }

  Tensor(const std::size_t _rank) : size(_rank, dimension), elements(getExtent())
  {
    std::fill(elements.begin(), elements.end(), value_type());
  }

  Tensor(const TensorSize& _size) : size(_size), elements(getExtent())
  {
    assert(size.getDimension() == dimension);
    std::fill(elements.begin(), elements.end(), value_type());
  }

  Tensor(const Tensor& t) : size(t.size), elements(t.elements)
  {
  }

  iterator begin()
  {
    return elements.begin();
  }

  const_iterator begin() const
  {
    return elements.begin();
  }

  iterator end()
  {
    return elements.end();
  }

  iterator end() const
  {
    return elements.end();
  }

  std::size_t getRank() const
  {
    return size.getRank();
  }

  std::size_t getDimension() const
  {
    return dimension;
  }
  
  TensorSize getSize() const
  {
    return size;
  }

  std::size_t getExtent() const
  {
    return size.getExtent();
  }
  
  value_type& operator()(const size_type index)
  {
    assert(getRank() == 1);
    assert(index < dimension);
    return elements[index];
  }

  const value_type operator()(const size_type index) const
  {
    assert(getRank() == 1);
    assert(index < dimension);
    return elements[index];
  }

  value_type& operator[](const size_type* const indices)
  {
    TensorIndex index(size);

    assert(getRank() == 0 || indices != NULL);

    for(std::size_t i=0; i<getRank(); ++i)
      index[i] = indices[i];

    return (*this)[index];
  }

  value_type& operator[](const TensorIndex& indices)
  {
    assert(indices.getSize() == size);
    return elements[TensorIndex::flatten(indices, row_major_tag())];
  }

  const value_type& operator[](const TensorIndex& indices) const
  {
    assert(indices.getSize() == size);
    return elements[TensorIndex::flatten(indices, row_major_tag())];
  }


  Tensor& operator*=(const value_type s)
  {
    std::transform(elements.begin(), elements.end(), elements.begin(), std::bind2nd(std::multiplies<value_type>(), s));
    return *this;
  }

  Tensor& operator/=(const value_type s)
  {
    std::transform(elements.begin(), elements.end(), elements.begin(), std::bind2nd(std::divides<value_type>(), s));
    return *this;
  }

  Tensor operator/(const value_type& t) const
  {
    Tensor result(*this);
    result/=t;
    return result;
  }

  Tensor operator*(const value_type s) const
  {
    Tensor x(*this);
    x*=s;
    return x;
  }

  Tensor& operator+=(const Tensor& t)
  {
    assert(size == t.size);
    std::transform(elements.begin(), elements.end(), t.elements.begin(), elements.begin(), std::plus<value_type>());
    return *this;
  }

  Tensor operator+(const Tensor& t) const
  {
    Tensor result(*this);
    result+=t;
    return result;
  }

  Tensor operator*(const Tensor& t) const
  {
    return outer_product(t);
  }

  Tensor outer_product(const Tensor& t) const
  {
    Tensor result(TensorSize(getRank() + t.getRank(), dimension));

    for(unsigned index=0; index<getExtent(); ++index)
      for(unsigned index2=0; index2<t.getExtent(); ++index2)
        result.elements[index * t.getExtent() + index2] = elements[index] * t.elements[index2];

    return result;
  }

  Tensor inner_product(const Tensor& t) const
  {
    assert(getRank() + t.getRank() >= 2);

    const std::size_t rank = getRank();
    const std::size_t iterationSize = pow(dimension, getRank()-1);
    const std::size_t tIterationSize = pow(dimension, t.getRank()-1);
    Tensor result(TensorSize(getRank() + t.getRank() - 2, dimension));

    for(unsigned index=0; index < iterationSize; ++index)
    {
      for(unsigned index2=0; index2 < tIterationSize; ++index2)
      {
        value_type sum = 0.0;
        for(unsigned sumIndex=0; sumIndex<dimension; ++sumIndex)
          sum += elements[index*rank + sumIndex] * t.elements[sumIndex*tIterationSize + index2];

        result.elements[index*tIterationSize + index2] = sum;
      }
    }

    return result;
  }

  Tensor colon_product(const Tensor& t) const
  {
    assert(getRank() == t.getRank());
    const TensorSize resultSize(0, dimension);
    Tensor result(resultSize);
    result.elements[0] = std::inner_product(elements.begin(), elements.end(), t.elements.begin(), value_type(0));
    return result;
  }

  void setElement(const std::size_t i, const Tensor& t)
  {
    assert(getRank() >= 1);
    assert(t.getRank() == getRank() - 1);
    assert(i < getDimension());

    const std::size_t elementSize = t.getExtent();
    for(unsigned index = 0; index < elementSize; ++index)
    {
      elements[i*elementSize + index] = t.elements[index]; 
    }
  }

  operator value_type()
  {
    assert(getRank() == 0 && "Attempt to convert non-rank 0 tensor to scalar");
    return elements[0];
  }

  void write(std::ostream& o) const
  {
    o << "Tensor: rank=" << size.getRank() << ", dimension=" << size.getDimension() << std::endl;

    const std::size_t extent = size.getExtent();
    for(std::size_t i=0; i<extent; ++i)
    {
      const TensorIndex index = TensorIndex::unflatten(size, i, row_major_tag());
      o << index << " = " << (*this)[index] << std::endl;
    }

    o << std::endl;
  }
};

template<std::size_t D, typename T> 
Tensor<D, T>  operator*(const typename Tensor<D, T>::value_type s, const Tensor<D, T>& t)
{
  return t*s;
}

template<std::size_t D, typename T> 
std::ostream& operator<<(std::ostream& o, const Tensor<D, T>& t)
{
  t.write(o);
  return o;
}

}

#endif
