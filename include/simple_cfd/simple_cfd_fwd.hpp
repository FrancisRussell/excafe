#ifndef SIMPLE_CFD_SIMPLE_CFD_FWD_HPP
#define SIMPLE_CFD_SIMPLE_CFD_FWD_HPP

#include <cstddef>

namespace cfd
{

// ID types
typedef std::size_t vertex_id;
typedef std::size_t cell_id;

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
template<typename C> class mesh_builder;
template<unsigned D> class mesh_geometry;
template<unsigned D> class mesh_geometry_impl;
template<unsigned int D> class vertex;
template<typename C> class finite_element;
template<typename T> class MeshFunction;
class TriangularCell;
class MeshTopology;
class MeshConnectivity;
class MeshEntity;
class MeshEntityIteratorGlobal;
class MeshEntityIteratorLocal;
class Polygon;

// Numeric types
class SparsityPattern;
class PETScMatrix;
class PETScVector;
class PETScKrylovSolver;
template<unsigned D, unsigned R, unsigned K, typename T> class TensorView;
template<unsigned int D, unsigned int R, typename T> class Tensor;
template<typename C> class FEMatrix;
template<typename C> class FEVector;

// Polynomial types
class Polynomial;
class Monomial;
class OptimisedPolynomial;
}

#endif
