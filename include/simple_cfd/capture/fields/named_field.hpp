#ifndef SIMPLE_CFD_CAPTURE_FIELDS_NAMED_FIELD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_NAMED_FIELD_HPP

#include <string>
#include "field.hpp"
#include "function_space.hpp"
#include "discrete_field_persistent.hpp"

namespace cfd
{

class NamedField
{
private:
  Field field;

public:
  NamedField()
  {
  }

  NamedField(const std::string& name, const FunctionSpace& functionSpace) : 
    field(new detail::DiscreteFieldPersistent(name, functionSpace))
  {
  }

  operator Field() const
  {
    return field;
  }
};

}

#endif
