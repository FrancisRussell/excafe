#ifndef SIMPLE_CFD_CAPTURE_SCENARIO_HPP
#define SIMPLE_CFD_CAPTURE_SCENARIO_HPP

#include <cstddef>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/mesh_function.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include "fields/element.hpp"
#include "fields/function_space.hpp"
#include "fields/function_space_mesh_function.hpp"

namespace cfd
{

template<std::size_t D>
class Scenario
{
private:
  static const std::size_t dimension = D;
  Mesh<dimension>& mesh;
  boost::ptr_vector< FiniteElement<dimension> > elements;

public:
  Scenario(Mesh<dimension>& _mesh) : mesh(_mesh)
  {
  }

  Element addElement(FiniteElement<dimension>* const e)
  {
    elements.push_back(e);
    return Element(elements.size() - 1);
  }
  
  FunctionSpace defineFunctionSpace(const Element& element, const Mesh<dimension>& _mesh)
  {
    assert(&mesh == &_mesh);
    return FunctionSpace(new detail::FunctionSpaceMeshFunction(element, MeshFunction<bool>(dimension, true)));
  }
};

}

#endif
