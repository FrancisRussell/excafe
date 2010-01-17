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
template<unsigned int D> class Function;
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
template<std::size_t> class Tensor;
template<std::size_t D> class DiscreteOperator;
template<std::size_t D> class DiscreteField;

// Polynomial types
class Polynomial;
template<typename V> class Monomial;
class OptimisedPolynomial;

// Boundary Conditions
template<std::size_t D> class BoundaryCondition3;
template<std::size_t D> class BoundaryConditionList;
template<std::size_t D> class BoundaryConditionTrivial;
}

#endif
