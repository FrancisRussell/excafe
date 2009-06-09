#ifndef SIMPLE_CFD_GENERIC_CELL_HPP
#define SIMPLE_CFD_GENERIC_CELL_HPP

#include<set>
#include<cstddef>
#include<simple_cfd_fwd.hpp>

namespace cfd
{

class GeneralCell
{
public:
  virtual std::set< std::set<std::size_t> > 
    getIncidentVertices(MeshTopology& topology, const MeshEntity& cellEntity, std::size_t d) const = 0;
  virtual std::size_t getDimension() const = 0;
};

}

#endif
