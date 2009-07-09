#ifndef SIMPLE_CFD_GENERIC_CELL_HPP
#define SIMPLE_CFD_GENERIC_CELL_HPP

#include<set>
#include<map>
#include<vector>
#include<cstddef>
#include<simple_cfd_fwd.hpp>

namespace cfd
{

class GeneralCell
{
public:
  virtual std::size_t getDimension() const = 0;
  virtual std::size_t getVerticesPerCell() const = 0;
  virtual std::size_t getLocalIndex(MeshTopology& topology, const std::size_t cid, const MeshEntity& entity) const = 0;
  virtual vertex<2> getLocalVertex(const std::size_t index) const = 0;
  virtual std::vector< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const = 0;
  virtual QuadraturePoints<2> getQuadrature(const std::size_t degree) const = 0;
  virtual double getJacobian(const mesh<TriangularCell>& m, const MeshEntity& entity, const vertex<2>& v) const = 0;
  virtual vertex<2> reference_to_physical(const mesh<TriangularCell>& m, const std::size_t cid, const vertex<2>& vertex) const = 0;
  virtual Tensor<2, 1, double> getFacetNormal(const mesh<TriangularCell>& m, const std::size_t cid, const
    std::size_t fid, const vertex<2>& v) const = 0;
  virtual ~GeneralCell() {}
};

}

#endif
