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
  virtual std::vector< std::set<std::size_t> > getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const = 0;
  virtual std::size_t getLocalIndex(MeshTopology& topology, const MeshEntity& entity, const std::size_t cid) const = 0;
  virtual std::size_t getDimension() const = 0;
  virtual std::size_t getVerticesPerCell() const = 0;
  virtual std::map<vertex<2>, double> getQuadrature(const mesh<TriangularCell>& m, const MeshEntity& entity) const = 0;
  virtual double getJacobian(const mesh<TriangularCell>& m, const MeshEntity& entity, const vertex<2>& v) const = 0;
  virtual Tensor<2, 1, double> getFacetNormal(const mesh<TriangularCell>& m, const std::size_t cid, const
    std::size_t fid, const vertex<2>& v) const = 0;
  virtual bool contains(const mesh<TriangularCell>& m, const std::size_t cid, const vertex<2>& v) const = 0;
  virtual ~GeneralCell() {}
};

}

#endif
