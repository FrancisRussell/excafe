#ifndef SIMPLE_CFD_MP_FLOAT_HPP
#define SIMPLE_CFD_MP_FLOAT_HPP

#include "mp_fwd.hpp"
#include <gmp.h>
#include <iosfwd>
#include <boost/operators.hpp>

namespace excafe
{

namespace mp
{

class Float : boost::totally_ordered<Float>,
              boost::arithmetic<Float>
{
private:
  friend class Integer;
  friend class Rational;

  mpf_t value;

public:
  Float();
  Float(const char* str);
  Float(double f);
  Float(const Float& f);
  Float(const Integer& i);
  Float(const Rational& r);
  Float& operator=(const Float& f);
  bool operator==(const Float& f) const;
  bool operator<(const Float& f) const;
  Float& operator+=(const Float& f);
  Float& operator-=(const Float& f);
  Float& operator*=(const Float& f);
  Float& operator/=(const Float& f);
  Float operator-() const;
  Float abs() const;
  Float pow(long exponent) const;
  void write(std::ostream& out) const;
  float toFloat() const;
  double toDouble() const;
  void swap(Float& f);
  ~Float();
};

std::ostream& operator<<(std::ostream& o, const Float& f);
Float abs(const Float& f);
Float pow(const Float& f, const int exponent);
Float pow(const Float& f, const long exponent);

}

}

#endif
