#ifndef SIMPLE_CFD_SIMPLE_CFD_FWD_HPP
#define SIMPLE_CFD_SIMPLE_CFD_FWD_HPP

#include <cstddef>

namespace cfd
{

// ID types
typedef std::size_t vertex_id;
typedef std::size_t cell_id;

// Mesh related types
template<std::size_t D> class Mesh;
template<unsigned D> class MeshGeometry;
template<unsigned D> class MeshGeometryImpl;
template<unsigned int D> class vertex;
template<std::size_t D> class FiniteElement;
template<unsigned int D> class SubDomain;
template<unsigned int D, unsigned int R> class Function;
template<typename T> class MeshFunction;
class MeshCell;
template<std::size_t> class GeneralCell;
class TriangularCell;
class MeshTopology;
class MeshConnectivity;
class MeshEntity;
class MeshEntityIteratorGlobal;
class MeshEntityIteratorLocal;
class Polygon;
template<std::size_t D> class CellVertices;

// Degree-of-freedom related types
template<std::size_t D> class Dof;
class DofAssociation;

// Quadrature related
class Quadrature;
template<std::size_t D> class QuadraturePoints;

// Numeric types
class SparsityPattern;
class PETScMatrix;
class PETScVector;
class PETScKrylovSolver;
template<unsigned D, unsigned R, unsigned K> class TensorView;
template<unsigned int D, unsigned int R> class Tensor;
template<typename C> class FEMatrix;
template<typename C> class FEVector;

// Polynomial types
class Polynomial;
class Monomial;
class OptimisedPolynomial;
}

#endif
