#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_TENSOR_ARRAY_REF_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_TENSOR_ARRAY_REF_HPP

#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/assert.hpp>
#include "index.hpp"
#include "tensor_fwd.hpp"
#include "tensor_array.hpp"

namespace cfd
{

namespace detail
{

class TensorArrayRef
{
private:
  typedef TensorArray element_t;
  boost::shared_ptr<element_t> value;

public:
  template<typename R>
  static TensorArrayRef cloneFrom(const R& _value)
  {
    BOOST_STATIC_ASSERT((boost::is_base_of<TensorArray, R>::value));
    return TensorArrayRef(new R(_value));
  }

  TensorArrayRef(element_t* const _value) : value(_value)
  {
  }

  TensorArrayRef(const TensorArrayRef& h) : value(h.value)
  {
  }

  element_t& operator*() const
  {
    return *value;
  }

  element_t* operator->() const
  {
    return value.get();
  }

  element_t* get() const
  {
    return value.get();
  }

  TensorArrayRef& operator=(const TensorArrayRef& rhs)
  {
    value = rhs.value;
    return *this;
  }

  bool operator==(const TensorArrayRef& rhs) const
  {
    return value == rhs.value;
  }

  bool operator<(const TensorArrayRef& rhs) const
  {
    return value < rhs.value;
  }

  ScalarPlaceholder operator[](const TensorIndex::constant_t i) const;
};

}

}

#endif
