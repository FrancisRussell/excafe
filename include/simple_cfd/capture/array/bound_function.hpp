#ifndef SIMPLE_CFD_CAPTURE_ARRAY_BOUND_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_BOUND_FUNCTION_HPP

#include <boost/variant.hpp>
#include "parameter_identifiers.hpp"

namespace cfd
{

namespace detail
{

class BoundFunction
{
private:
  std::map< ArrayIndexID, boost::variant<std::size_t, ArrayIndexID> > arrayIndexBindings;
  std::map< TensorIndexID, boost::variant<std::size_t, TensorIndexID> > tensorIndexBindings;
  std::map< ScalarID, BoundFunction> scalarBindings;

public:
  BoundFunction()
  {
  }
};

}

}

#endif
