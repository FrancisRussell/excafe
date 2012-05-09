#ifndef EXCAFE_UTIL_LAZY_COPY_HPP
#define EXCAFE_UTIL_LAZY_COPY_HPP

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace excafe
{

namespace util
{

/*
   This class supplies copy-on-write semantics for an arbitrary class. The class needs to be well-behaved
   with respect to const-correctness otherwise disaster may ensue, especially in the multi-threaded case.
   In addition, although the ref() and cref() methods may look more appealing than using than operator*(),
   avoid them where possible. The finer granularity of ref() and cref() make a mistake like passing an
   iterator from a transparently cloned object to its clone more likely.
*/

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

  LazyCopy(const LazyCopy& l) : value(l.value)
  {
  }

  template<typename T1>
  explicit LazyCopy(const T1& v) : value(boost::make_shared<value_type>(v))
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
  void swap(excafe::util::LazyCopy<T>& a, excafe::util::LazyCopy<T>& b)
  {
    a.swap(b);
  }
}

#endif
