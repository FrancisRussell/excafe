#ifndef SIMPLE_CFD_UTIL_LAZY_COPY_HPP
#define SIMPLE_CFD_UTIL_LAZY_COPY_HPP

#include <boost/shared_ptr.hpp>

namespace cfd
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

  struct ValueHolder
  {
    ValueHolder() {}
    ValueHolder(const value_type& _value) : value(_value) {}

    value_type value;
  };

  boost::shared_ptr<ValueHolder> holder;

public:
  LazyCopy() : holder(new ValueHolder())
  {
  }

  LazyCopy(const value_type& v) : holder(new ValueHolder(v))
  {
  }

  LazyCopy(const LazyCopy& l) : holder(l.holder)
  {
  }

  LazyCopy& operator=(const LazyCopy& l)
  {
    holder = l.holder;
    return *this;
  }

  const value_type& cref() const
  {
    return holder->value;
  }

  value_type& ref()
  {
    /* Is this problematic in the multi-threaded case? */
    if (!holder.unique())
      holder = boost::shared_ptr<ValueHolder>(new ValueHolder(holder->value));

    return holder->value;
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
    holder.swap(l.holder);
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
