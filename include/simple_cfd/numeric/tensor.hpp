#ifndef SIMPLE_CFD_NUMERIC_TENSOR_HPP
#define SIMPLE_CFD_NUMERIC_TENSOR_HPP

#include <simple_cfd_fwd.hpp>
#include <boost/array.hpp>
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

template<unsigned D, unsigned R, unsigned K, typename T>
class TensorView
{
private:
  typedef T value_type;
  typedef std::size_t size_type;

  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  static const unsigned int knownIndicesCount = K;

  Tensor<D, R, T>& tensor;
  boost::array<std::size_t, knownIndicesCount>  knownIndices;

public:
  TensorView(Tensor<D, R, T>& t, size_type* indices) : tensor(t)
  {
    std::copy(indices, indices+knownIndicesCount, knownIndices);
  }

  TensorView<D, R, K+1, T> operator()(const size_type index)
  {
     std::size_t newKnownIndices[knownIndicesCount+1];
     std::copy(knownIndices, knownIndices+knownIndicesCount, newKnownIndices);
     newKnownIndices[knownIndicesCount] = index;

     return TensorView<D, R, K+1, T>(tensor, newKnownIndices);
  }
};

template<unsigned D, unsigned R, typename T>
class TensorView<D, R, R, T>
{
private:
  typedef T value_type;
  typedef std::size_t size_type;
  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  
  Tensor<D, R, T>& tensor;
  std::size_t knownIndices[rank];

public:
  TensorView(Tensor<D, R, T>& t, const size_type* const indices) : tensor(t)
  {
    std::copy(indices, indices+rank, knownIndices);
  }

  T toScalar() const
  {
    return tensor.toScalar();
  }
};


template<unsigned int D, unsigned int R, typename T>
class Tensor
{
public:
  typedef T value_type;
  typedef std::size_t size_type;
  static const unsigned int dimension = D;
  static const unsigned int rank = R;
  static const size_type size = detail::Power<dimension, rank>::value;

  template<unsigned int D_, unsigned int R_, typename T_>
  friend class Tensor;

private:
  boost::array<value_type, size>  elements;

public:
  Tensor()
  {
    std::fill(elements.begin(), elements.end(), T());
  }

  Tensor(const Tensor& t) : elements(t.elements)
  {
  }

  TensorView<D, R, 1, T> operator()(const size_type index)
  {
    return TensorView<D, R, 1, T>(*this, &index);
  }

  T& operator[](const size_type* const indices)
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

  Tensor& operator*=(const T s)
  {
    std::transform(elements.begin(), elements.end(), elements.begin(), std::bind2nd(std::multiplies<T>(), s));
    return *this;
  }

  Tensor operator*(const T s) const
  {
    Tensor x(*this);
    x*=s;
    return x;
  }

  template<unsigned int R2>
  Tensor<D, R+R2, T> operator*(const Tensor<D, R2, T>& t) const
  {
    Tensor<D, R+R2, T> result;

    for(unsigned index=0; index<detail::Power<D, R>::value; ++index)
      for(unsigned index2=0; index2<detail::Power<D, R2>::value; ++index2)
        result.elements[index * detail::Power<D, R2>::value + index2] = elements[index] * t.elements[index2];

    return result;
  }

  Tensor<D, 0, T> colon_product(const Tensor& t) const
  {
    Tensor<D, 0, T> result;
    result.elements[0] = std::inner_product(elements.begin(), elements.end(), t.elements.begin(), T(0));
    return result;
  }

  T toScalar() const
  {
    assert(rank == 0);
    return elements[0];
  }
};

template<unsigned int D, unsigned int R, typename T> 
Tensor<D, R, T>  operator*(const T s, const Tensor<D, R, T>& t)
{
  return t*s;
}

}

#endif
