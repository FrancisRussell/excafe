#ifndef SIMPLE_CFD_CAPTURE_NUMERIC_CLN_WRAPPER_HPP
#define SIMPLE_CFD_CAPTURE_NUMERIC_CLN_WRAPPER_HPP

#include <ostream>
#include <boost/operators.hpp>
#include <cln/cln.h>

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

  void swap(CLNWrapper& b)
  {
    std::swap(value, b.value);
  }

  void write(std::ostream& o) const
  {
    o << value;
  }

  operator double() const
  {
    return cln::double_approx(value);
  }

  operator cln::cl_F() const
  {
    return value;
  }
};

template<std::size_t P>
std::ostream& operator<<(std::ostream& o, const CLNWrapper<P> c)
{
  c.write(o);
  return o;
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
