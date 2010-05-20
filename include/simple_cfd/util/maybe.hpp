#ifndef SIMPLE_CFD_UTIL_MAYBE_HPP
#define SIMPLE_CFD_UTIL_MAYBE_HPP

#include <ostream>
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace util
{

class Nothing {};

template<typename T>
class Maybe
{
public:
  typedef T value_type;

  bool nothing;
  value_type val;

public:
  Maybe() : nothing(true)
  {
  }

  Maybe(const Nothing& n) : nothing(true)
  {
  }

  Maybe(const Maybe& m) : nothing(m.nothing), val(m.val)
  {
  }

  Maybe(const value_type& v) : nothing(false), val(v)
  {
  }

  bool hasValue() const
  {
    return !nothing;
  }

  bool isNothing() const
  {
    return nothing;
  }

  const value_type value() const
  {
    if (nothing)
      CFD_EXCEPTION("Tried to retrieve value from a Maybe with nothing inside.");
    else
      return val;
  }

  void write(std::ostream& o) const
  {
    if (nothing)
      o << "Nothing";
    else
      o << val;
  }
};

}

}

namespace std
{

template<typename T>
std::ostream& operator<<(std::ostream& o, const cfd::util::Maybe<T>& m)
{
  m.write(o);
  return o;
}

}

#endif
