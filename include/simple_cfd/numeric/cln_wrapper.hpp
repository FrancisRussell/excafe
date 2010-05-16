#ifndef SIMPLE_CFD_CAPTURE_NUMERIC_CLN_WRAPPER_HPP
#define SIMPLE_CFD_CAPTURE_NUMERIC_CLN_WRAPPER_HPP

#include <ostream>
#include <boost/operators.hpp>
#include <cln/cln.h>
#include "cast.hpp"

namespace cfd
{

template<std::size_t P>
class CLNWrapper : boost::totally_ordered<CLNWrapper<P>,
                   boost::arithmetic<CLNWrapper<P>
                   > >
{
public:
  static const std::size_t precision = P;

private:
  cln::cl_F value;

  template<typename T>
  static cln::cl_F convert(const T& s)
  {
    static const cln::float_format_t format = cln::float_format(precision);
    return cln::cl_float(s, format);
  }

public:
  CLNWrapper() : value(convert(0.0))
  {
  }

  CLNWrapper(const double s) : value(convert(s))
  {
  }
  
  CLNWrapper(const float s) : value(convert(s))
  {
  }

  // We need to call convert on CLN values because they may have a different precision
  CLNWrapper(const cln::cl_F& c) : value(convert(c))
  {
  }

  CLNWrapper(const CLNWrapper& c) : value(c.value)
  {
  }

  bool operator<(const CLNWrapper& c) const
  {
    return value < c.value;
  }

  bool operator==(const CLNWrapper& c) const
  {
    return value == c.value;
  }

  CLNWrapper& operator+=(const CLNWrapper& c)
  {
    value += c.value;
    return *this;
  }

  CLNWrapper& operator-=(const CLNWrapper& c)
  {
    value -= c.value;
    return *this;
  }

  CLNWrapper& operator*=(const CLNWrapper& c)
  {
    value *= c.value;
    return *this;
  }

  CLNWrapper& operator/=(const CLNWrapper& c)
  {
    value /= c.value;
    return *this;
  }

  double toDouble() const
  {
    return cln::double_approx(value);
  }

  float toFloat() const
  {
    return cln::float_approx(value);
  }

  void swap(CLNWrapper& b)
  {
    std::swap(value, b.value);
  }

  CLNWrapper pow(const int y) const
  {
    return convert(cln::expt(value, y));
  }

  CLNWrapper floor() const
  {
    return convert(cln::floor1(value));
  }

  CLNWrapper ceil() const
  {
    return convert(cln::ceiling1(value));
  }

  void write(std::ostream& o) const
  {
    o << value;
  }
};

template<std::size_t P>
CLNWrapper<P> pow(const CLNWrapper<P>& x, const CLNWrapper<P>& y)
{
  return x.pow(y);
}

template<std::size_t P>
CLNWrapper<P> pow(const CLNWrapper<P>& x, const int y)
{
  return x.pow(y);
}

template<std::size_t P>
CLNWrapper<P> floor(const CLNWrapper<P>& x)
{
  return x.floor();
}

template<std::size_t P>
CLNWrapper<P> ceil(const CLNWrapper<P>& x)
{
  return x.ceil();
}

template<std::size_t P>
std::ostream& operator<<(std::ostream& o, const CLNWrapper<P> c)
{
  c.write(o);
  return o;
}

namespace detail
{

template<std::size_t P>
struct RawConverter<CLNWrapper<P>, float>
{
  static float low_level_convert(const CLNWrapper<P>& s) 
  { 
    return s.toFloat();
  }
};


template<std::size_t P>
struct RawConverter<CLNWrapper<P>, double>
{
  static double low_level_convert(const CLNWrapper<P>& s) 
  { 
    return s.toDouble();
  }
};

}

}

namespace std
{

template<std::size_t P>
void swap(cfd::CLNWrapper<P>& a, cfd::CLNWrapper<P>& b)
{
  a.swap(b);
}

}

#endif
