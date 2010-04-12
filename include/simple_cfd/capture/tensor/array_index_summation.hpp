#ifndef SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_INDEX_SUMMATION_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_ARRAY_INDEX_SUMMATION_HPP

#include <set>
#include "tensor_array.hpp"
#include "tensor_size.hpp"
#include "tensor_array_ref.hpp"
#include <simple_cfd/exception.hpp>
#include <boost/foreach.hpp>

namespace cfd
{

namespace detail
{

class ArrayIndexSummation : public TensorArray
{
private:
  TensorArrayRef ref;
  std::set<ArrayIndexVariable> variables;

  void validate() const
  {
    if (variables.empty()) CFD_EXCEPTION("ArrayIndexSummation require at least one index.");
    const std::size_t limit = variables.begin()->getLimit();

    BOOST_FOREACH(const ArrayIndexVariable& var, variables)
    {
      if (limit != var.getLimit()) CFD_EXCEPTION("ArrayIndexSummation requires all variables have same limit.");
    }
  }

public:
  ArrayIndexSummation(const TensorArrayRef& _ref, const ArrayIndexVariable& v) : ref(_ref)
  {
    variables.insert(v);
    validate();
  }

  ArrayIndexSummation(const TensorArrayRef& _ref, const std::set<ArrayIndexVariable>& _variables) : 
    ref(_ref), variables(_variables)
  {
    validate();
  }

  virtual TensorSize getTensorSize() const
  {
    return ref->getTensorSize();
  }

  virtual TensorArrayRef derivative(const ScalarPlaceholder& x) const 
  {
    ArrayIndexSummation result(*this);
    result.ref = result.ref->derivative(x);

    return TensorArrayRef::cloneFrom(result);
  }
};

}

}

#endif
