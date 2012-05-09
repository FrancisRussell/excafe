#ifndef EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_MESH_HPP
#define EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_MESH_HPP

#include "element.hpp"
#include "function_space_visitor.hpp"

namespace excafe
{

namespace detail
{

class FunctionSpaceMesh : public FunctionSpaceExpr
{
private:
  Element element;  

public:
  FunctionSpaceMesh(const Element& _element) :
    element(_element)
  {
  }

  void accept(FunctionSpaceVisitor& v)
  {
    v.visit(*this);
  }

  Element getElement() const
  {
    return element;
  }
};

}

}

#endif
