#ifndef SIMPLE_CFD_CAPTURE_SCENARIO_HPP
#define SIMPLE_CFD_CAPTURE_SCENARIO_HPP

#include <cstddef>
#include <string>
#include <map>
#include <set>
#include <cassert>
#include <utility>
#include <ostream>
#include <fstream>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/foreach.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/mesh_function.hpp>
#include <simple_cfd/dof_map.hpp>
#include <simple_cfd/discrete_field.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/boundary_condition_list.hpp>
#include <simple_cfd/local_assembly_matrix.hpp>
#include "solve_operation.hpp"
#include "dimensionless_scenario.hpp"
#include "fields/element.hpp"
#include "fields/function_space.hpp"
#include "fields/function_space_expr.hpp"
#include "fields/function_space_mesh.hpp"
#include "fields/named_field.hpp"
#include "fields/boundary_condition.hpp"
#include "evaluation/evaluation_fwd.hpp"
#include "evaluation/function_space_resolver.hpp"
#include "evaluation/boundary_condition_builder.hpp"
#include "evaluation/assembly_optimising_visitor.hpp"

namespace cfd
{

template<std::size_t D>
class Scenario : public detail::DimensionlessScenario
{
private:
  static const std::size_t dimension = D;
  typedef vertex<dimension> vertex_type;
  typedef typename DofMap<dimension>::dof_t dof_t;
  typedef detail::FunctionSpaceExpr* function_space_ptr;
  typedef detail::LocalAssemblyMatrix<dimension, detail::ScalarPlaceholder::optimised_expression_t> local_assembly_matrix_t;

  Mesh<dimension>& mesh;
  boost::ptr_vector< FiniteElement<dimension> > elements;
  std::map< function_space_ptr, DofMap<dimension> > functionSpaceMap;
  std::map< std::string, DiscreteField<dimension> > persistentFields;
  std::map< const detail::OperatorAssembly*, local_assembly_matrix_t > optimisedCellIntegrals;
  std::vector< DiscreteField<dimension> > boundaryValues;

  boost::tuple<double, double, double> getValue(const std::size_t cid, const vertex_type& vertex, 
    const DiscreteField<dimension>& field) const
  {
    assert(dimension <= 3);
    assert(field.getRank() <= 1);

    boost::tuple<double, double, double> value(0.0, 0.0, 0.0);

    for(unsigned dof=0; dof<field.getElement()->spaceDimension(); ++dof)
    {
      Tensor<dimension> basis = field.getElement()->evaluateTensor(dof, vertex);
      const dof_t fieldDof = dof_t(field.getElement(), cid, dof);

      double valueCoeff;
      field.getValues(1u, &fieldDof, &valueCoeff);

      if (field.getRank() == 0)
      {
        boost::get<0>(value) += basis * valueCoeff;
      }
      else
      {
        if (dimension >= 1)
          boost::get<0>(value) += basis(0) * valueCoeff;
        if (dimension >= 2)
         boost::get<1>(value) += basis(1) * valueCoeff;
        if (dimension >= 3)
         boost::get<2>(value) += basis(2) * valueCoeff;
      }
    }
    return value;
  }


  void renderFields(std::ostream& out) const
  {
    typedef typename std::map< std::string, DiscreteField<dimension> >::const_iterator named_field_map_iter;
    static const std::size_t vtk_dimension = 3;    
    assert(dimension <= 3);

    std::vector< vertex<dimension> > vertices;
    std::map< std::string, std::vector< boost::tuple<double, double, double> > > values;

    for(typename Mesh<dimension>::global_iterator vIter(mesh.global_begin(0)); vIter!=mesh.global_end(0); ++vIter)
    {
      const vertex<dimension> v(mesh.getVertex(vIter->getIndex()));
      vertices.push_back(v);

      const std::size_t cid = mesh.getContainingCell(*vIter);
      const MeshEntity localVertexEntity = mesh.getLocalEntity(cid, *vIter);
      const vertex<dimension> localVertex = mesh.getLocalCoordinate(cid, localVertexEntity.getIndex());

      for (named_field_map_iter fieldIter(persistentFields.begin()); fieldIter != persistentFields.end(); ++fieldIter)
        values[fieldIter->first].push_back(getValue(cid, localVertex, fieldIter->second));
    }

    for(typename Mesh<dimension>::global_iterator eIter(mesh.global_begin(1)); eIter!=mesh.global_end(1); ++eIter)
    {
      const std::size_t cid = mesh.getContainingCell(*eIter);
      const std::vector<std::size_t> vertexIndices(mesh.getIndices(*eIter, 0));
      assert(vertexIndices.size() == 2);
      
      const MeshEntity v1Entity = mesh.getLocalEntity(cid, MeshEntity(0, vertexIndices[0]));
      const MeshEntity v2Entity = mesh.getLocalEntity(cid, MeshEntity(0, vertexIndices[1]));

      const vertex<dimension> localVertex((mesh.getLocalCoordinate(cid, v1Entity.getIndex()) + 
        mesh.getLocalCoordinate(cid, v2Entity.getIndex()))/2.0);

      vertices.push_back(mesh.referenceToPhysical(cid, localVertex));

      for (named_field_map_iter fieldIter(persistentFields.begin()); fieldIter != persistentFields.end(); ++fieldIter)
        values[fieldIter->first].push_back(getValue(cid, localVertex, fieldIter->second));
    }

    out << "# vtk DataFile Version 2.0" << std::endl;
    out << "Simple Navier-Stokes Solver" << std::endl;
    out << "ASCII" << std::endl;
    out << "DATASET POLYDATA" << std::endl;
    out << "POINTS " << vertices.size() << " DOUBLE " << std::endl;

    for(std::size_t point = 0; point < vertices.size(); ++point)
    {
      const vertex<dimension> v(vertices[point]);

      for(std::size_t index=0; index < vtk_dimension; ++index)
      {
        if (index < dimension)
          out << v[index];
        else
          out << "0";

        out << " ";
      }
      out << std::endl;
    }

    out << "POLYGONS ";
    out << mesh.numEntities(dimension) << " ";
    out << mesh.numRelations(dimension, 0) + mesh.numEntities(dimension) << std::endl; 

    for(typename Mesh<dimension>::global_iterator cIter(mesh.global_begin(dimension)); cIter!=mesh.global_end(dimension); ++cIter)
    {
      const std::vector<std::size_t> vIndices(mesh.getIndices(*cIter, 0));
      out << vIndices.size();

      for(std::size_t v=0; v<vIndices.size(); ++v)
      {
        out << " " << vIndices[v];
      }
      out << std::endl;
    }

    out << "POINT_DATA " << vertices.size() << std::endl;
    for (typename std::map< std::string, DiscreteField<dimension> >::const_iterator fieldIter(persistentFields.begin());
      fieldIter != persistentFields.end(); ++fieldIter)
    {
      const bool isScalarField = (fieldIter->second.getRank() == 0);

      out << (isScalarField ? "SCALARS" : "VECTORS") << " " << fieldIter->first << " DOUBLE" << std::endl;
      
      if (isScalarField)
        out << "LOOKUP_TABLE default" << std::endl;

      for(std::size_t point = 0; point < values[fieldIter->first].size(); ++point)
      {
        out << boost::get<0>(values[fieldIter->first][point]) << " ";
        if (!isScalarField)
        {
          out << boost::get<1>(values[fieldIter->first][point]) << " ";
          out << boost::get<2>(values[fieldIter->first][point]);      
        }
        out << std::endl;
      }
    }
  }

  void resolveFunctionSpace(detail::FunctionSpaceExpr& f)
  {
    // TODO: make me efficent
    std::set<function_space_ptr> functionSpaceSet;
    functionSpaceSet.insert(&f);
    resolveFunctionSpaces(functionSpaceSet);
  }

public:
  Scenario(Mesh<dimension>& _mesh) : mesh(_mesh)
  {
  }

  Mesh<dimension>& getMesh()
  {
    return mesh;
  }

  const Mesh<dimension>& getMesh() const
  {
    return mesh;
  }

  Element addElement(FiniteElement<dimension>* const e)
  {
    elements.push_back(e);
    return Element(elements.size() - 1);
  }

  const FiniteElement<dimension>& getElement(const Element& e) const
  {
    assert(e.getIndex() < elements.size());
    return elements[e.getIndex()];
  }
  
  FunctionSpace defineFunctionSpace(const Element& element, const Mesh<dimension>& _mesh)
  {
    assert(&mesh == &_mesh);
    return FunctionSpace(new detail::FunctionSpaceMesh(element));
  }

  NamedField defineNamedField(const std::string& name, const FunctionSpace functionSpace)
  {
    if (persistentFields.find(name) != persistentFields.end())
    {
      CFD_EXCEPTION("Attempted to create two named fields with the same name.");
    }

    resolveFunctionSpace(*functionSpace.getExpr());
    NamedField field(name, functionSpace);
    persistentFields.insert(std::make_pair(name, DiscreteField<dimension>(getDofMap(*functionSpace.getExpr()))));
    return field;
  }

  virtual void resolveFunctionSpaces(const std::set<function_space_ptr> functionSpaces)
  {
    detail::FunctionSpaceResolver<dimension> functionSpaceResolver(*this, functionSpaceMap);

    for(std::set<function_space_ptr>::const_iterator fsIter(functionSpaces.begin()); fsIter!=functionSpaces.end(); ++fsIter)
      (*fsIter)->accept(functionSpaceResolver);
  }

  virtual void optimiseLocalAssemblies(const detail::DiscreteExprSet<detail::discrete_operator_tag>& operators)
  {
    detail::AssemblyOptimisingVisitor<dimension> visitor(*this);

    BOOST_FOREACH(detail::OperatorExpr& expr, std::make_pair(operators.begin_expr(), operators.end_expr()))
    {
      expr.accept(visitor);
    }
  }

  void setOptimisedCellIntegral(const detail::OperatorAssembly& assembly, const local_assembly_matrix_t& matrix)
  {
    optimisedCellIntegrals.insert(std::make_pair(&assembly, matrix));
  }

  std::map<MeshEntity, local_assembly_matrix_t> getOptimisedCellIntegral(const detail::OperatorAssembly& assembly) const
  {
    const MeshEntity cellEntity(dimension, 0);
    typename std::map<const detail::OperatorAssembly*, local_assembly_matrix_t>::const_iterator matIter =
      optimisedCellIntegrals.find(&assembly);
    
    if (matIter != optimisedCellIntegrals.end())
    {
      std::map<MeshEntity, local_assembly_matrix_t> matrices;
      matrices.insert(std::make_pair(cellEntity, matIter->second));
      return matrices;
    }
    else
    {
      CFD_EXCEPTION("Unable to locate optimised cell integral for specified assembly operator.");
    }
  }

  const DofMap<dimension>& getDofMap(detail::FunctionSpaceExpr& e) const
  {
    const typename std::map< function_space_ptr, DofMap<dimension> >::const_iterator mapIter = functionSpaceMap.find(&e);
    assert(mapIter != functionSpaceMap.end());
    return mapIter->second;
  }

  DiscreteField<dimension> getNamedValue(const detail::DiscreteFieldPersistent& p)
  {
    const typename std::map< std::string, DiscreteField<dimension> >::iterator fieldIter = persistentFields.find(p.getName());
    assert(fieldIter != persistentFields.end());
    return fieldIter->second;
  }

  DiscreteField<dimension> getBoundaryField(const BoundaryCondition& b)
  {
    return boundaryValues[b.getIndex()];
  }

  void setNamedValue(const std::string& name, const DiscreteField<dimension>& v)
  {
    DiscreteField<dimension> newValue(v);
    const typename std::map< std::string, DiscreteField<dimension> >::iterator fieldIter = persistentFields.find(name);
    assert(fieldIter != persistentFields.end());
    std::swap(fieldIter->second, newValue);
  }

  BoundaryCondition addBoundaryCondition(const FunctionSpace& f, const BoundaryConditionList<dimension>& condition)
  {
    resolveFunctionSpace(*f.getExpr());

    detail::BoundaryConditionBuilder<dimension> builder(mesh);
    boundaryValues.push_back(builder.getBoundaryValues(getDofMap(*f.getExpr()), condition));
    return BoundaryCondition(boundaryValues.size() - 1);
  }

  SolveOperation newSolveOperation()
  {
    return SolveOperation(*this);
  }

  void outputFieldsToFile(const std::string& filename) const
  {
    std::ofstream outFile(filename.c_str());
    renderFields(outFile);
    outFile.close();
  }

  void execute(SolveOperation& o)
  {
    o.executeDimensionTemplated<dimension>(*this);
  }
};

}

#endif
