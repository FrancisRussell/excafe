#ifndef SIMPLE_CFD_CAPTURE_FIELDS_NAMED_FIELD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_NAMED_FIELD_HPP

#include "field.hpp"
#include "function_space.hpp"
#include "discrete_field_persistent.hpp"
#include <string>
#include <cassert>

namespace cfd
{

class NamedField
{
private:
  bool defined;
  std::string name;
  Field field;

public:
  NamedField() : defined(false)
  {
  }

  NamedField(const std::string& _name, const FunctionSpace& functionSpace) : 
    defined(true), name(_name), field(new detail::DiscreteFieldPersistent(name, functionSpace.getExpr()))
  {
    assert(!name.empty());
  }

  bool operator<(const NamedField& f) const
  {
    return name < f.name;
  }

  bool operator==(const NamedField& f) const
  {
    return name == f.name;
  }

  Field getField()
  {
    return field;
  }

  operator Field() const
  {
    return field;
  }
};

}

#endif
