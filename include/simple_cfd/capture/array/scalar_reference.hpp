#ifndef SIMPLE_CFD_CAPTURE_ARRAY_SCALAR_REFERENCE_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_SCALAR_REFERENCE_HPP

#include <boost/variant.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/operators.hpp>
#include "array_traits.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include "free_tensor_array.hpp"
#include "tensor_function_constant.hpp"

namespace cfd
{

namespace detail
{

class ScalarReference : boost::equality_comparable<ScalarReference>
{
public:
  typedef FreeTensorArray unbound_ref_t;
  typedef TensorFunctionConstant::value_ptr bound_ref_t;

private:
  typedef boost::variant<unbound_ref_t, bound_ref_t> ref_t;

  ArrayIndex<param_tag> arrayIndex;
  TensorIndex<param_tag> tensorIndex;
  ref_t ref;

public:
  ScalarReference(const ArrayIndex<param_tag>& _arrayIndex, const TensorIndex<param_tag>& _tensorIndex, 
    const unbound_ref_t& _ref) :
    arrayIndex(_arrayIndex), tensorIndex(_tensorIndex), ref(_ref)
  {
  }

  ScalarReference(const ArrayIndex<param_tag>& _arrayIndex, const TensorIndex<param_tag>& _tensorIndex, 
    const ScalarReference& _ref) :
    arrayIndex(_arrayIndex), tensorIndex(_tensorIndex), ref(_ref.ref)
  {
  }

  bool operator==(const ScalarReference& s) const
  {
    return arrayIndex == s.arrayIndex && 
           tensorIndex == s.tensorIndex &&
           ref == s.ref;
  }

  bool operator<(const ScalarReference& s) const
  {
    return boost::make_tuple(arrayIndex, tensorIndex, ref) < boost::make_tuple(arrayIndex, tensorIndex, ref);
  }

  bool isBound() const
  {
    return boost::get<bound_ref_t>(&ref) != NULL;
  }

  unbound_ref_t getFreeTensorArray() const
  {
    assert(!isBound());
    return boost::get<unbound_ref_t>(ref);
  }

  ArrayIndex<param_tag> getArrayIndex() const
  {
    return arrayIndex;
  }

  TensorIndex<param_tag> getTensorIndex() const
  {
    return tensorIndex;
  }

  bool isParameterised() const
  {
    return arrayIndex.isParameterised() || tensorIndex.isParameterised(); 
  }
};

}

}
#endif
