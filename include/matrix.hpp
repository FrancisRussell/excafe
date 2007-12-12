#ifndef SIMPLE_CFD_MATRIX_HPP
#define SIMPLE_CFD_MATRIX_HPP

#include<cassert>
#include<map>
#include<utility>
#include<ostream>
#include<iomanip>

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
    assert(width >= 0);
    assert(height >= 0);
  }

  matrix() : width(0), height(0)
  {
  }

  value_type& operator()(const size_type x, const size_type y)
  {
    assert(x>=0 && x<width);
    assert(y>=0 && y<height);
    return values[std::make_pair(x, y)];
  }

  const value_type operator()(const size_type x, const size_type y) const
  {
    assert(x>=0 && x<width);
    assert(y>=0 && y<height);

    const typename std::map<std::pair<size_type, size_type>, value_type>::const_iterator valueIter(values.find(std::make_pair(x, y)));
    if (valueIter == values.end())
      return 0;
    else
      return valueIter->second;
  }

  size_type getWidth() const
  {
    return width;
  }

  size_type getHeight() const
  {
    return height;
  }

  void resize(const unsigned w, const unsigned h)
  {
    width = w;
    height = h;
    
    assert(width >= 0);
    assert(height >= 0);
  }
};

}

template<typename T>
std::ostream& operator<<(std::ostream& o, const cfd::matrix<T>& m)
{
  const unsigned rows = m.getWidth();
  const unsigned cols = m.getHeight();

  o << "Matrix: " << rows << " * " << cols << std::endl << "[" << std::endl;
  for(unsigned row=0; row<rows; ++row)
  {
    o << "[";
    for(unsigned col=0; col<cols; ++col)
    {
      o << std::setw(9) << std::fixed;
      o << m(row,col);
      if (col!=cols-1)
        o <<", ";
    }
    o << "]" << std::endl;
  }
  o << "]";
  return o << std::endl;
}

#endif
