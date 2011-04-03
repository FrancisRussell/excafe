#ifndef SIMPLE_CFD_UTIL_LAZY_COPY_HPP
#define SIMPLE_CFD_UTIL_LAZY_COPY_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace cfd
{

namespace util
{

template<typename T>
class LazyCopy
{
public:
  typedef T value_type;

private:
  bool operator==(const LazyCopy&) const;
  boost::shared_ptr<value_type> value;

public:
  LazyCopy() : value(boost::make_shared<value_type>())
  {
  }

  explicit LazyCopy(const value_type& v) : value(boost::make_shared<value_type>(v))
  {
  }

  LazyCopy(const LazyCopy& l) : value(l.value)
  {
  }

  LazyCopy& operator=(const LazyCopy& l)
  {
    value = l.value;
    return *this;
  }

  const value_type& cref() const
  {
    return *value;
  }

  value_type& ref()
  {
    /* Is this problematic in the multi-threaded case? */
    if (!value.unique())
      value = boost::make_shared<value_type>(*value);

    return *value;
  }

  const value_type& operator*() const
  {
    return cref();
  }

  value_type& operator*()
  {
    return ref();
  }

  const value_type* operator->() const
  {
    return &cref();
  }

  value_type* operator->()
  {
    return &ref();
  }

  void swap(LazyCopy& l)
  {
    value.swap(l.value);
  }
};

}

}

namespace std
{
  template<typename T>
  void swap(cfd::util::LazyCopy<T>& a, cfd::util::LazyCopy<T>& b)
  {
    a.swap(b);
  }
}

#endif
