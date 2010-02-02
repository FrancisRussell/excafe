#ifndef SIMPLE_CFD_CAPTURE_ARRAY_SCALAR_REFERENCE_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_SCALAR_REFERENCE_HPP

#include <boost/variant.hpp>
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "free_tensor_array.hpp"
#include "tensor_function_constant_holder.hpp"

namespace cfd
{

namespace detail
{

class ScalarReference
{
public:
  typedef FreeTensorArray unbound_ref_t;
  typedef TensorFunctionConstant::function_ptr bound_ref_t;

private:
  typedef boost::variant<unbound_ref_t, bound_ref_t> ref_t;

  ArrayIndex arrayIndex;
  TensorIndex tensorIndex;
  ref_t ref;

public:
  ScalarReference(const ArrayIndex& _arrayIndex, const TensorIndex& _tensorIndex, const unbound_ref_t& _ref) :
    arrayIndex(_arrayIndex), tensorIndex(_tensorIndex), _ref(ref)
  {
  }
};

}

}
#endif
