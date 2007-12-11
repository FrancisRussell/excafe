#ifndef SIMPLE_CFD_MATRIX_HPP
#define SIMPLE_CFD_MATRIX_HPP

#include<cassert>
#include<map>
#include<utility>

namespace cfd
{

template<typename T>
class matrix
{
public:
  typedef T value_type;
  typedef int size_type;
private:
  size_type width;
  size_type height;
  std::map<std::pair<size_type, size_type>, value_type> values;

public:
  matrix(const size_type w, const size_type h) : width(w), height(h)
  {
    assert(width>=0);
    assert(height>=0);
  }

  value_type& operator()(const size_type x, const size_type y)
  {
    assert(x>=0 && x<width);
    assert(y>=0 && y<height);
    return values[std::pair<size_type, size_type>(x, y)];
  }

  const value_type& operator()(const size_type x, const size_type y) const
  {
    assert(x>=0 && x<width);
    assert(y>=0 && y<height);
    return values[std::pair<size_type, size_type>(x, y)];
  }

};

}

#endif
