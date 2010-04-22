#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_CELL_VERTICES_PLACEHOLDER_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_CELL_VERTICES_PLACEHOLDER_HPP

#include "scalar_placeholder.hpp"
#include "cell_vertex_component.hpp"

namespace cfd
{

namespace detail
{

template<std::size_t D>
class CellVerticesPlaceholder
{
public:
  static const std::size_t dimension = D;
  typedef Tensor<dimension, assembly_polynomial_t> tensor_t;

  tensor_t operator[](const std::size_t i) const
  {
    const TensorSize resultSize(1, dimension);
    tensor_t result(resultSize);
    TensorIndex index(resultSize);

    for(std::size_t d=0; d<dimension; ++d)
    {
      index[0] = d;
      result[index] = ScalarPlaceholder(CellVertexComponent(i, d));
    }

    return result;
  }
};

}

}

#endif
