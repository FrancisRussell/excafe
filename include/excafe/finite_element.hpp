#ifndef EXCAFE_FINITE_ELEMENT_HPP
#define EXCAFE_FINITE_ELEMENT_HPP

#include <vector>
#include <set>
#include <cstddef>
#include <map>
#include <ostream>
#include "excafe_fwd.hpp"
#include "cell_manager.hpp"
#include "numeric/tensor.hpp"
#include "numeric/polynomial.hpp"
#include "capture/assembly/assembly_fwd.hpp"
#include "capture/assembly/position_placeholder.hpp"
#include "capture/assembly/scalar_placeholder.hpp"

namespace excafe
{

template<std::size_t D>
class FiniteElement
{
public:
  static const std::size_t dimension = D;
  typedef vertex<dimension> vertex_type;
  typedef Dof<dimension> dof_t;
  typedef typename CellManager::ref<dimension>::general cell_ref_t;
  typedef detail::ScalarPlaceholder::expression_t expression_t;
  typedef Tensor<dimension, expression_t> tensor_expr_t;

  virtual std::size_t getRank() const = 0;
  virtual std::size_t getDimension() const = 0;
  virtual unsigned getTensorIndex(const Mesh<dimension>& mesh, const std::size_t cid, const std::size_t dof) const = 0;
  virtual unsigned spaceDimension() const = 0; // Number of basis functions
  virtual std::vector< std::set<dof_t> > resolveIdenticalDofs(const Mesh<dimension>& m, 
    const MeshEntity& entity, 
    const std::set<dof_t>& dofsOnEntity) const = 0;
  virtual vertex_type getDofCoordinateLocal(const std::size_t dof) const = 0;
  virtual vertex_type getDofCoordinateGlobal(const Mesh<dimension>& m, const cell_id cid, const std::size_t dof) const = 0;
  virtual std::set< Dof<dimension> > getDofsOnEntity(MeshTopology& topology, const cell_id cid, const MeshEntity& entity) const = 0;
  virtual tensor_expr_t getBasis(const std::size_t i, const detail::PositionPlaceholder& pos) const = 0;
  virtual Tensor<dimension> evaluateTensor(const std::size_t i, const vertex_type& vRef) const = 0;
  virtual Tensor<dimension> evaluateDivergence(const std::size_t i, const vertex_type& vRef) const = 0;
  virtual Tensor<dimension> evaluateGradient(const std::size_t i, const vertex_type& vRef) const = 0;
  virtual cell_ref_t getCell() const = 0;
  virtual void write(std::ostream& o) const = 0;
  virtual ~FiniteElement() {}
};

}

namespace std
{

template<std::size_t D>
std::ostream& operator<<(std::ostream& o, const excafe::FiniteElement<D>& f) 
{
  f.write(o);
  return o;
}

}

#endif
