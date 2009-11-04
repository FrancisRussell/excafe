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

template<std::size_t D>
class Tensor
{
public:
  typedef double value_type;
  typedef std::size_t size_type;
  static const std::size_t  dimension = D;
  std::size_t rank;
  size_type size;

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
  Tensor() : rank(0), size(pow(dimension, rank)), elements(size)
  {
    std::fill(elements.begin(), elements.end(), value_type());
  }

  Tensor(const std::size_t _rank) : rank(_rank), size(pow(dimension, rank)), elements(size)
  {
    std::fill(elements.begin(), elements.end(), value_type());
  }

  Tensor(const Tensor& t) : rank(t.rank), size(t.size), elements(t.elements)
  {
  }

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  std::size_t numValues() const
  {
    return size;
  }
  
  value_type& operator()(const size_type index)
  {
    assert(rank == 1);
    assert(index < dimension);
    return elements[index];
  }

  value_type& operator[](const size_type* const indices)
  {
    assert(rank == 0 || indices != NULL);

    std::size_t index = 0;
  
    for(unsigned i = 0; i<rank; ++i)
    {
      assert(indices[i] < dimension);
      index *= dimension;
      index += indices[i];
    }

    assert(index < size);
    return elements[index];
  }

  Tensor& operator*=(const value_type s)
  {
    std::transform(elements.begin(), elements.end(), elements.begin(), std::bind2nd(std::multiplies<value_type>(), s));
    return *this;
  }

  Tensor operator*(const value_type s) const
  {
    Tensor x(*this);
    x*=s;
    return x;
  }

  Tensor& operator+=(const Tensor& t)
  {
    assert(rank == t.rank);
    std::transform(elements.begin(), elements.end(), t.elements.begin(), elements.begin(), std::plus<double>());
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
    Tensor result(rank + t.rank);

    for(unsigned index=0; index<size; ++index)
      for(unsigned index2=0; index2<t.size; ++index2)
        result.elements[index * t.size + index2] = elements[index] * t.elements[index2];

    return result;
  }

  Tensor inner_product(const Tensor& t) const
  {
    assert(rank + t.rank >= 2);

    const std::size_t iterationSize = pow(dimension, rank-1);
    const std::size_t tIterationSize = pow(dimension, t.rank-1);
    Tensor result(rank + t.rank - 2);

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

  value_type colon_product(const Tensor& t) const
  {
    assert(rank == t.rank);
    return std::inner_product(elements.begin(), elements.end(), t.elements.begin(), value_type(0));
  }

  operator value_type()
  {
    assert(rank == 0 && "Attempt to convert non-rank 0 tensor to scalar");
    return elements[0];
  }
};

template<unsigned int D> 
Tensor<D>  operator*(const typename Tensor<D>::value_type s, const Tensor<D>& t)
{
  return t*s;
}

}

#endif
