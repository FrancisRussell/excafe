#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_MESH_FUNCTION_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_MESH_FUNCTION_HPP

#include <simple_cfd/mesh_function.hpp>
#include "element.hpp"
#include "function_space_visitor.hpp"

namespace cfd
{

namespace detail
{

class FunctionSpaceMeshFunction : public FunctionSpaceExpr
{
private:
  Element element;  
  MeshFunction<bool> subdomain;

public:
  FunctionSpaceMeshFunction(const Element& _element, const MeshFunction<bool>& _subdomain) :
    element(_element), subdomain(_subdomain)
  {
  }

  void accept(FunctionSpaceVisitor& v)
  {
    v.visit(*this);
  }
};

}

}

#endif
