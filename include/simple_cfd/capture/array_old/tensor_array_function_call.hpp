#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_CALL_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_ARRAY_FUNCTION_CALL_HPP

#include <boost/variant.hpp>
#include <boost/tuple/tuple.hpp>
#include "tensor_function.hpp"
#include "array_index.hpp"
#include "tensor_index.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

class TensorArrayFunctionCall
{
private:
  class undefined_tag {};
  typedef boost::tuple<ArrayIndex<param_tag>, TensorIndex<param_tag>, TensorFunction::ref> call_tuple_t;

  boost::variant<undefined_tag, call_tuple_t> call;

  const call_tuple_t& getTuple() const
  {
    try
    {
      return boost::get<const call_tuple_t&>(call);
    }
    catch(boost::bad_get& e)
    {
      CFD_EXCEPTION("Tried to access values of an uninitialised tensor function call.");
    }
  }

  call_tuple_t& getTuple()
  {
    try
    {
      return boost::get<call_tuple_t&>(call);
    }
    catch(boost::bad_get& e)
    {
      CFD_EXCEPTION("Tried to access values of an uninitialised tensor function call.");
    }
  }

public:
  TensorArrayFunctionCall() : call(undefined_tag())
  {
  }

  TensorArrayFunctionCall(const ArrayIndex<param_tag>& _arrayIndex, const TensorIndex<param_tag>& _tensorIndex, 
    const TensorFunction::ref _function) : call(call_tuple_t(_arrayIndex, _tensorIndex, _function))
  {
  }

  ArrayIndex<param_tag> getArrayIndex() const
  {
    return boost::get<0>(getTuple());
  }

  TensorIndex<param_tag> getTensorIndex() const
  {
    return boost::get<1>(getTuple());
  }

  TensorFunction::ref getFunction() const
  {
    return boost::get<2>(getTuple());
  }

  void setArrayIndex(const ArrayIndex<param_tag>& _arrayIndex)
  {
    boost::get<0>(getTuple()) = _arrayIndex;
  }

  void setTensorIndex(const TensorIndex<param_tag>& _tensorIndex)
  {
    boost::get<1>(getTuple()) = _tensorIndex;
  }

  void setFunction(const TensorFunction::ref _function)
  {
    boost::get<2>(getTuple()) = _function;
  }
};

}

}

#endif
