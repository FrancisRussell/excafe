#ifndef SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_CONSTANT_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_TENSOR_FUNCTION_CONSTANT_HPP

#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class TensorFunctionConstant
{
public:
  typedef boost::shared_ptr<TensorFunctionConstant> value_ptr;
};

}

}
#endif
