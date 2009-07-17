#ifndef SIMPLE_CFD_NUMERIC_TENSOR_HPP
#define SIMPLE_CFD_NUMERIC_TENSOR_HPP

#include <simple_cfd_fwd.hpp>
#include <boost/array.hpp>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <cassert>
#include <algorithm>
#include <functional>
#include <numeric>

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

template<unsigned D, unsigned R, unsigned K>
class TensorView
{
private:
  typedef typename Tensor<D, R>::value_type value_type;
  typedef std::size_t size_type;

  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  static const unsigned int knownIndicesCount = K;

  Tensor<D, R>& tensor;
  boost::array<std::size_t, knownIndicesCount>  knownIndices;

public:
  TensorView(Tensor<D, R>& t, const size_type* const indices) : tensor(t)
  {
    std::copy(indices, indices+knownIndicesCount, knownIndices.begin());
  }

  TensorView<D, R, K+1> operator()(const size_type index)
  {
     std::size_t newKnownIndices[knownIndicesCount+1];
     std::copy(knownIndices.begin(), knownIndices.end(), newKnownIndices);
     newKnownIndices[knownIndicesCount] = index;

     return TensorView<D, R, K+1>(tensor, newKnownIndices);
  }
};

template<unsigned D, unsigned R>
class TensorView<D, R, R>
{
private:
  typedef typename Tensor<D, R>::value_type value_type;
  typedef std::size_t size_type;
  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  
  Tensor<D, R>& tensor;
  std::size_t knownIndices[rank];

public:
  TensorView(Tensor<D, R>& t, const size_type* const indices) : tensor(t)
  {
    std::copy(indices, indices+rank, knownIndices);
  }

  value_type toScalar() const
  {
    return tensor[knownIndices];
  }

  value_type& operator=(const value_type value)
  {
    return tensor[knownIndices] = value;
  }
};


template<unsigned int D, unsigned int R>
class Tensor
{
public:
  typedef double value_type;
  typedef std::size_t size_type;
  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  static const size_type size = detail::Power<dimension, rank>::value;

  template<unsigned int D_, unsigned int R_>
  friend class Tensor;

private:
  std::vector<value_type>  elements;

public:
  Tensor() : elements(size)
  {
    std::fill(elements.begin(), elements.end(), value_type());
  }

  Tensor(const Tensor& t) : elements(t.elements)
  {
  }

  TensorView<D, R, 1> operator()(const size_type index)
  {
    return TensorView<D, R, 1>(*this, &index);
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

  template<unsigned int R2>
  Tensor<D, R+R2> operator*(const Tensor<D, R2>& t) const
  {
    Tensor<D, R+R2> result;

    for(unsigned index=0; index<detail::Power<D, R>::value; ++index)
      for(unsigned index2=0; index2<detail::Power<D, R2>::value; ++index2)
        result.elements[index * detail::Power<D, R2>::value + index2] = elements[index] * t.elements[index2];

    return result;
  }

  template<unsigned int R2>
  Tensor<D, R+R2-2> inner_product(const Tensor<D, R2>& t) const
  {
    Tensor<D, R+R2-2> result;

    for(unsigned index=0; index<detail::Power<D, R-1>::value; ++index)
    {
      for(unsigned index2=0; index2<detail::Power<D, R2-1>::value; ++index2)
      {
        value_type sum = 0.0;
        for(unsigned sumIndex=0; sumIndex<D; ++sumIndex)
          sum += elements[index*rank + sumIndex] * t.elements[sumIndex*detail::Power<D, R2-1>::value + index2];

        result.elements[index*detail::Power<D, R2-1>::value + index2] = sum;
      }
    }

    return result;
  }

  Tensor<D, 0> colon_product(const Tensor& t) const
  {
    Tensor<D, 0> result;
    result.elements[0] = std::inner_product(elements.begin(), elements.end(), t.elements.begin(), value_type(0));
    return result;
  }

  value_type toScalar() const
  {
    assert(rank == 0);
    return elements[0];
  }
};

template<unsigned int D, unsigned int R> 
Tensor<D, R>  operator*(const typename Tensor<D, R>::value_type s, const Tensor<D, R>& t)
{
  return t*s;
}

}

#endif
