#ifndef EXCAFE_NUMERIC_SMALL_VECTOR_HPP
#define EXCAFE_NUMERIC_SMALL_VECTOR_HPP

#include <cstddef>
#include <algorithm>
#include <functional>
#include <boost/bind.hpp>
#include <boost/array.hpp>
#include <boost/operators.hpp>
#include <excafe/exception.hpp>
#include "numeric_fwd.hpp"

namespace excafe
{

template<std::size_t NumRows, typename T>
class SmallVector : boost::additive<SmallVector<NumRows, T>,
                    boost::multiplicative<SmallVector<NumRows, T>, T
                    > >
{
public:
  typedef std::size_t size_type;
  static const std::size_t N = NumRows;
  typedef T value_type;
  typedef typename boost::array<value_type, N>::iterator iterator;
  typedef typename boost::array<value_type, N>::const_iterator const_iterator;

private:
  boost::array<value_type, N> values;

public:
  SmallVector()
  {
  }

  SmallVector(const SmallVector& v) : values(v.values)
  {
  }

  SmallVector(const Tensor<N, value_type>& v)
  {
    if (v.getRank() != 1)
    {
      CFD_EXCEPTION("Tried to construct a vector from a non-rank 1 tensor.");
    }
    else
    {
      std::copy(v.begin(), v.end(), begin());
    }
  }

  SmallVector& operator=(const SmallVector& v)
  {
    values = v.values;
    return *this;
  }

  iterator begin() 
  {
    return values.begin();
  }

  iterator end()
  {
    return values.end();
  }

  const_iterator begin() const
  {
    return values.begin();
  }

  const_iterator end() const
  {
    return values.end();
  }

  size_type numRows() const
  {
    return N;
  }

  SmallVector& operator*=(const value_type s)
  {
    std::transform(begin(), end(), begin(), boost::bind(std::multiplies<value_type>(), _1, s));
    return *this;
  }

  SmallVector& operator/=(const value_type s)
  {
    std::transform(begin(), end(), begin(), boost::bind(std::divides<value_type>(), _1, s));
    return *this;
  }

  SmallVector& operator+=(const SmallVector& s)
  {
    std::transform(begin(), end(), s.begin(), begin(), std::plus<value_type>());
    return *this;
  }

  SmallVector& operator-=(const SmallVector& s)
  {
    std::transform(begin(), end(), s.begin(), begin(), std::minus<value_type>());
    return *this;
  }

  value_type& operator[](const size_type i)
  {
    return values[i];
  }

  const value_type& operator[](const size_type i) const
  {
    return values[i];
  }
};

}

#endif
