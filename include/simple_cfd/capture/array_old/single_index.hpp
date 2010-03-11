#ifndef SIMPLE_CFD_CAPTURE_ARRAY_SINGLE_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_SINGLE_INDEX_HPP

#include <cstddef>
#include <boost/variant.hpp>
#include "parameter_identifiers.hpp"

namespace cfd
{

namespace detail
{

template<typename P>
class SingleIndex
{
public:
  typedef std::size_t constant_t;
  typedef P           parameter_t;

private:
  boost::variant<constant_t, parameter_t> value;

public:
  SingleIndex() : value(0)
  {
  }

  SingleIndex(const constant_t c) : value(c)
  {
  }

  SingleIndex(const parameter_t& p) : value(p)
  {
  }

  bool isConstant() const
  {
    return boost::get<constant_t>(&value) != NULL;
  }

  bool isParameter() const
  {
    return boost::get<parameter_t>(&value) != NULL;
  }

  constant_t getConstant() const
  {
    return boost::get<constant_t>(value);
  }

  parameter_t getParameter() const
  {
    return boost::get<parameter_t>(value);
  }

  bool operator==(const SingleIndex& s) const
  {
    return value == s.value;
  }

  bool operator<(const SingleIndex& s) const
  {
    return value < s.value;
  }
};

}

}

#endif
