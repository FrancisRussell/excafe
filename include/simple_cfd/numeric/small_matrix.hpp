#ifndef SIMPLE_CFD_NUMERIC_SMALL_MATRIX_HPP
#define SIMPLE_CFD_NUMERIC_SMALL_MATRIX_HPP

#include <cstddef>
#include <algorithm>
#include <numeric>
#include <functional>
#include <boost/array.hpp>
#include <boost/bind.hpp>
#include <boost/operators.hpp>
#include "numeric_fwd.hpp"

namespace cfd
{

template<std::size_t NumRows, std::size_t NumCols, typename T>
class SmallMatrix : boost::additive<SmallMatrix<NumRows, NumCols, T>,
                    boost::multiplicative<SmallMatrix<NumRows, NumCols, T>, T
                    > >
{
public:
  typedef std::size_t size_type;
  static const std::size_t M = NumRows;
  static const std::size_t N = NumCols;
  typedef T value_type;
  typedef SmallVector<N, value_type> OneD;
  typedef typename boost::array<OneD, M>::iterator iterator;
  typedef typename boost::array<OneD, M>::const_iterator const_iterator;

private:
  boost::array<OneD, M> values;

public:
  SmallMatrix()
  {
  }

  SmallMatrix(const SmallMatrix& m) : values(m.values)
  {
  }

  value_type& operator()(const size_type row, const size_type col)
  {
    return values[row][col];
  }

  const value_type& operator()(const size_type row, const size_type col) const
  {
    return values[row][col];
  }

  OneD& operator[](const size_type i)
  {
    return values[i];
  }

  const OneD& operator[](const size_type i) const
  {
    return values[i];
  }

  iterator begin() 
  {
    return values.begin();
  }

  const_iterator begin() const
  {
    return values.begin();
  }

  iterator end() 
  {
    return values.end();
  }

  const_iterator end() const
  {
    return values.end();
  }

  size_type numRows() const
  {
    return M;
  }

  size_type numCols() const
  {
    return N;
  }

  SmallMatrix& operator*=(const value_type s)
  {
    std::transform(begin(), end(), begin(), boost::bind(std::multiplies<value_type>(), _1, s));
    return *this;
  }

  SmallMatrix& operator/=(const value_type s)
  {
    std::transform(begin(), end(), begin(), boost::bind(std::divides<value_type>(), _1, s));
    return *this;
  }

  SmallMatrix& operator+=(const SmallMatrix& s)
  {
    std::transform(begin(), end(), s.begin(), begin(), std::plus<value_type>());
    return *this;
  }

  SmallMatrix& operator-=(const SmallMatrix& s)
  {
    std::transform(begin(), end(), s.begin(), begin(), std::minus<value_type>());
    return *this;
  }

  SmallVector<M, value_type> operator*(const SmallVector<N, value_type>& x) const
  {
    SmallVector<M, value_type> b;

    for(size_t r=0; r<M; ++r)
    {
      b[r] = std::inner_product((*this)[r].begin(), (*this)[r].end(), x.begin(), value_type());
    }

    return b;
  }
};
template<typename T>
T determinant(const SmallMatrix<1, 1, T>& m)
{
  return m(0,0);
}


template<typename T>
T determinant(const SmallMatrix<2, 2, T>& m)
{
  return m(0,0)*m(1,1) - m(0,1)*m(1,0);
}

template<typename T>
T determinant(const SmallMatrix<3, 3, T>& m)
{
  return  m(0,0)*m(1,1)*m(2,2) + m(0,1)*m(1,2)*m(2,0) + m(0,3)*m(1,0)*m(2,1)
         -m(0,0)*m(1,2)*m(2,1) - m(1,1)*m(1,0)*m(2,2) - m(2,2)*m(1,1)*m(2,0);
}


}

#endif
