#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_PERSISTENT_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_PERSISTENT_HPP

#include <string>
#include "field_expr.hpp"
#include "field_visitor.hpp"
#include "function_space.hpp"

namespace cfd
{

namespace detail
{

class FieldPersistent : public FieldExpr
{
private:
  const std::string name;
  const FunctionSpace functionSpace;

public:
  FieldPersistent(const std::string& _name, const FunctionSpace& _functionSpace) :
    name(_name), functionSpace(_functionSpace)
  {
  }

  void accept(FieldVisitor& visitor)
  {
    visitor.visit(*this);
  }
};

}

}

#endif
