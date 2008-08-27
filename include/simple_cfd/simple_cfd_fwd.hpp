#ifndef SIMPLE_CFD_SIMPLE_CFD_FWD_HPP
#define SIMPLE_CFD_SIMPLE_CFD_FWD_HPP

namespace cfd
{

// ID types
typedef int vertex_id;
typedef int cell_id;

// Cell shapes
enum shape
{
  triangle
};

// Traits call for shapes
template<shape s>
struct shape_dimensions
{
};

template<>
struct shape_dimensions<triangle>
{
  static const int dimension = 2;
};

// Mesh related types
template<typename C> class mesh;
template<unsigned D> class mesh_geometry;
template<unsigned D> class mesh_geometry_impl;
template<shape s> class cell;
template<unsigned int D> class vertex;
template<typename C> class finite_element;

// Basis function types
struct evaluated_basis
{
  double value;
  double dx;
  double dy;
};

// Numeric types
class SparsityPattern;
class PETScMatrix;
class PETScVector;
class PETScKrylovSolver;
template<unsigned D, unsigned R, unsigned K, typename T> class TensorView;
template<unsigned int D, unsigned int R, typename T> class Tensor;

}

#endif
