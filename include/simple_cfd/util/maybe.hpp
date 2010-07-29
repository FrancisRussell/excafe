#ifndef SIMPLE_CFD_UTIL_MAYBE_HPP
#define SIMPLE_CFD_UTIL_MAYBE_HPP

#include <ostream>
#include <algorithm>
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

  bool operator==(const Maybe& m) const
  {
    if (nothing == m.nothing)
    {
      return nothing || (val == m.val);
    }
    else
    {
      return false;
    }
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

  void swap(Maybe<T>& m)
  {
    std::swap(nothing, m.nothing);
    std::swap(val, m.val);
  }
};

}

}

namespace std
{

template<typename T>
void swap(cfd::util::Maybe<T>& a, cfd::util::Maybe<T>& b)
{
  a.swap(b);
}

}

template<typename T>
std::ostream& operator<<(std::ostream& o, const cfd::util::Maybe<T>& m)
{
  m.write(o);
  return o;
}

#endif
