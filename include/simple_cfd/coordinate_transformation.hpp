#ifndef SIMPLE_CFD_COORDINATE_TRANSFORMATION_HPP
#define SIMPLE_CFD_COORDINATE_TRANSFORMATION_HPP

#include <cstddef>
#include <cassert>
#include <utility>
#include <map>
#include "numeric/small_vector.hpp"
#include "numeric/small_matrix.hpp"
#include "numeric/cast.hpp"
#include <simple_cfd/capture/assembly/scalar_placeholder.hpp>
#include <simple_cfd/capture/assembly/position_placeholder.hpp>
#include <simple_cfd/capture/assembly/cell_vertices_placeholder.hpp>
#include <simple_cfd/capture/assembly/cell_vertex_component.hpp>

namespace cfd
{

template<std::size_t FromDimension, std::size_t ToDimension>
class TransformationBase
{
protected:
  typedef double value_type;
  typedef detail::ScalarPlaceholder::expression_t expression_t;
  static const std::size_t from_dim = FromDimension;
  static const std::size_t to_dim = ToDimension;

  SmallVector<to_dim, expression_t> transform;

  template<std::size_t D>
  static void addPositionPlaceholders(expression_t::value_map& map, const vertex<D>& v)
  {
    detail::PositionPlaceholder positionPlaceholder;
    for(std::size_t d=0; d<D; ++d)
    {
      map.bind(positionPlaceholder[d], v[d]);
    }
  }

  template<std::size_t D>
  static void addCellVertexPlaceholders(expression_t::value_map& map, const CellVertices<D>& vertices)
  {
    for(std::size_t vid=0; vid<vertices.size(); ++vid)
    {
      for(std::size_t d=0; d<from_dim; ++d)
      {
        map.bind(detail::ScalarPlaceholder(detail::CellVertexComponent(vid, d)), vertices[vid][d]);
      }
    }
  }

public:
  SmallVector<to_dim, expression_t> getTransformed() const
  {
    return this->transform;
  }

  SmallMatrix<to_dim, from_dim, expression_t> getJacobian() const
  {
    detail::PositionPlaceholder positionPlaceholder;
    SmallMatrix<to_dim, from_dim, expression_t> jacobian;

    for(std::size_t row=0; row<to_dim; ++row)
    {
      for(std::size_t col=0; col<from_dim; ++col)
      {
        jacobian(row, col) = this->transform[row].derivative(positionPlaceholder[col]);
      }
    }
    return jacobian;
  }

  expression_t getScalingFactor() const
  {
    return determinant(getJacobian());
  }

  SmallVector<to_dim, expression_t> getNormal() const
  {
    CFD_EXCEPTION("Invalid for this transformation.");
  }
};

template<std::size_t FromDimension, std::size_t ToDimension>
class GlobalTransformation : public TransformationBase<FromDimension, ToDimension>
{
private:
  typedef TransformationBase<FromDimension, ToDimension> base_t;
  typedef typename CellManager::ref<ToDimension>::general cell_ref_t;

public:
  typedef typename base_t::value_type value_type;
  typedef typename base_t::expression_t expression_t;
  static const std::size_t from_dim = base_t::from_dim;
  static const std::size_t to_dim = base_t::to_dim;

public:
  GlobalTransformation(const FiniteElement<to_dim>& element, const GeneralCell<to_dim>& cell)
  {
    detail::CellVerticesPlaceholder<to_dim> cellVertices;
    detail::PositionPlaceholder position;
    const std::size_t numVertices = cell.numEntities(0);
    assert(element.spaceDimension() == numVertices);

    for(std::size_t i=0; i<numVertices; ++i)
    {
      this->transform += element.getBasis(i, position).toScalar() * 
                      SmallVector<to_dim, expression_t>(cellVertices[i]);
    }
  }

  vertex<to_dim> getTransformed(const vertex<from_dim>& v, const CellVertices<to_dim>& vertices) const
  {
    typename expression_t::value_map placeholderMap;
    addPositionPlaceholders(placeholderMap, v);
    addCellVertexPlaceholders(placeholderMap, vertices);

    vertex<to_dim> transformed;

    for(std::size_t i=0; i<this->transform.numRows(); ++i)
      transformed[i] = cfd::numeric_cast<value_type>(this->transform[i].evaluate(placeholderMap));

    return transformed;
  }

  SmallMatrix<to_dim, from_dim, value_type> getJacobian(const vertex<from_dim>& v, const CellVertices<to_dim>& vertices) const
  {
    const SmallMatrix<to_dim, from_dim, expression_t> jacobianExpression = base_t::getJacobian();

    typename expression_t::value_map valueMap;
    addPositionPlaceholders(valueMap, v);
    addCellVertexPlaceholders(valueMap, vertices);

    SmallMatrix<to_dim, from_dim, value_type> jacobian;
    for(std::size_t row=0; row<to_dim; ++row)
    {
      for(std::size_t col=0; col<from_dim; ++col)
      {
        jacobian(row, col) = cfd::numeric_cast<value_type>(jacobianExpression(row, col).evaluate(valueMap));
      }
    }
    return jacobian;
  }

  value_type getScalingFactor(const vertex<from_dim>& v, const CellVertices<to_dim>& vertices) const
  {
    return determinant(getJacobian(v, vertices));
  }

  SmallVector<to_dim, double> getNormal(const vertex<from_dim>& v, const CellVertices<to_dim>& vertices) const
  {
    CFD_EXCEPTION("Invalid for this transformation.");
  }
};


}

#endif
