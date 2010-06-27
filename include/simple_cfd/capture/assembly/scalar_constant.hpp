#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_CONSTANT_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_SCALAR_CONSTANT_HPP

#include <cstddef>
#include <utility>
#include <ostream>
#include <cln/cln.h>
#include <simple_cfd/numeric/cast.hpp>

namespace cfd
{

namespace detail
{

class ScalarConstant
{
private:
  cln::cl_F constant;

public:
  ScalarConstant(const double s) : constant(s)
  {
  }

  bool operator==(const ScalarConstant& c) const
  {
    return constant == c.constant;
  }

  bool operator<(const ScalarConstant& c) const
  {
    return constant < c.constant;
  }

  void write(std::ostream& o) const
  {
    o << "const(" << constant << ")";
  }

  operator double() const
  {
    return cfd::numeric_cast<double>(constant);
  }

  operator float() const
  {
    return cfd::numeric_cast<float>(constant);
  }

  operator cln::cl_F() const
  {
    return constant;
  }
};

std::ostream& operator<<(std::ostream& o, const ScalarConstant& c);

}

}

#endif
