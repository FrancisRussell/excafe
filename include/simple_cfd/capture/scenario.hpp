#ifndef SIMPLE_CFD_CAPTURE_SCENARIO_HPP
#define SIMPLE_CFD_CAPTURE_SCENARIO_HPP

#include <cstddef>
#include <string>
#include <cassert>
#include <boost/ptr_container/ptr_vector.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/mesh_function.hpp>
#include <simple_cfd/dof_map.hpp>
#include "solve_operation.hpp"
#include "dimensionless_scenario.hpp"
#include "fields/element.hpp"
#include "fields/function_space.hpp"
#include "fields/function_space_expr.hpp"
#include "fields/function_space_mesh_function.hpp"
#include "fields/named_field.hpp"
#include "evaluation/function_space_resolver.hpp"

namespace cfd
{

template<std::size_t D>
class Scenario : public detail::DimensionlessScenario
{
private:
  static const std::size_t dimension = D;
  typedef detail::FunctionSpaceExpr* function_space_ptr;

  Mesh<dimension>& mesh;
  boost::ptr_vector< FiniteElement<dimension> > elements;
  std::map< function_space_ptr, DofMap<dimension> > functionSpaceMap;
  std::vector<NamedField> persistentFields;

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

  NamedField defineNamedField(const std::string& name, const FunctionSpace functionSpace)
  {
    NamedField field(name, functionSpace);
    persistentFields.push_back(field);
    return field;
  }

  virtual void resolveFunctionSpaces(const std::set<function_space_ptr> functionSpaces)
  {
    // TODO: implement me!
    detail::FunctionSpaceResolver<dimension> functionSpaceResolver(functionSpaceMap);
    assert(false);
  }

  SolveOperation newSolveOperation()
  {
    return SolveOperation(*this);
  }

  void outputFieldsToFile(const std::string& filename) const
  {
    // TODO: implement me!
    assert(false);
  }

  void execute(SolveOperation& o)
  {
    o.executeDimensionTemplated<dimension>();
  }
};

}

#endif
