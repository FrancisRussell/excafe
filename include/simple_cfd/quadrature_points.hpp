#ifndef SIMPLE_CFD_QUADRATURE_POINTS_HPP
#define SIMPLE_CFD_QUADRATURE_POINTS_HPP

#include <cstddef>
#include <vector>
#include <utility>
#include <map>
#include <cassert>
#include <simple_cfd_fwd.hpp>
#include <vertex.hpp>

namespace cfd
{

template<std::size_t D>
class QuadraturePoints
{
public:
  static const std::size_t dimension = D;

private:
  typedef std::vector< std::pair<vertex<dimension>, double> > qpoints_t;
  std::map<MeshEntity, qpoints_t> quadratureMap;

public:
  typedef typename qpoints_t::const_iterator const_iterator;
  typedef const_iterator iterator;

  QuadraturePoints() 
  {
  }

  void setQuadrature(const MeshEntity& e, const std::map<vertex<dimension>, double> weights)
  {
    assert(e.getDimension() <= dimension);
    // We only have one cell per cell, so if we have quadrature for MeshEntity(dimension, k)
    // for k>0, something is wrong
    assert(e.getDimension() < dimension || e.getIndex() == 0);

    quadratureMap[e] = qpoints_t(weights.begin(), weights.end());
  }

  bool hasQuadrature(const MeshEntity& e)
  {
    return quadratureMap.find(e) != quadratureMap.end();
  }

  const_iterator begin(const MeshEntity& e) const
  {
    assert(e.getDimension() <= dimension);
    assert(e.getDimension() < dimension || e.getIndex() == 0);

    const typename std::map<MeshEntity, qpoints_t>::const_iterator qIter = quadratureMap.find(e);
    assert(qIter != quadratureMap.end());
    return qIter->begin();
  }

  const_iterator end(const MeshEntity& e) const
  {
    assert(e.getDimension() <= dimension);
    assert(e.getDimension() < dimension || e.getIndex() == 0);

    const typename std::map<MeshEntity, qpoints_t>::const_iterator qIter = quadratureMap.find(e);
    assert(qIter != quadratureMap.end());
    return qIter->end();
  }
};

}

#endif
