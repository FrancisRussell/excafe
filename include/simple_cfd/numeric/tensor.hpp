#ifndef SIMPLE_CFD_NUMERIC_TENSOR_HPP
#define SIMPLE_CFD_NUMERIC_TENSOR_HPP

#include <simple_cfd_fwd.hpp>
#include <boost/array.hpp>
#include <algorithm>
#include <cstddef>
#include <cassert>

namespace cfd
{

namespace detail
{

template<unsigned  X, unsigned Y>
struct Power
{
  static const unsigned int value = Power<X*Y, Y-1>::value;
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
  std::size_t knownIndices[knownIndicesCount];

public:
  TensorView(Tensor<D, R, T>& t, size_type* indices) : tensor(t)
  {
    std::copy(indices, indices+knownIndicesCount, knownIndices);
  }

  TensorView<D, R, K+1, T> operator[](const size_type index)
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
  TensorView(Tensor<D, R, T>& t, size_type* indices) : tensor(t)
  {
    std::copy(indices, indices+rank, knownIndices);
  }

  T toScalar() const
  {
    tensor.toScalar();
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
  static const size_type size = detail::Power<D, R>::value;

private:
  boost::array<value_type, size>  elements;

public:
  Tensor()
  {
    std::fill(elements.begin(), elements.end(), T());
  }

  TensorView<D, R, 0, T> operator[](const size_type index)
  {
    return TensorView<D, R, 0, T>(this, NULL);
  }

  T& operator[](const size_type* const indices)
  {
    std::size_t index = 0;
  
    for(int i = 0; i<rank; ++i)
    {
      assert(indices[i] < dimension);
      index *= dimension;
      index += indices[i];
    }

    return elements[index];
  }

  T toScalar() const
  {
    assert(rank == 0);
    return elements[0];
  }
};

}

#endif
