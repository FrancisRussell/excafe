#ifndef SIMPLE_CFD_FINITE_ELEMENT_HPP
#define SIMPLE_CFD_FINITE_ELEMENT_HPP

#include "simple_cfd_fwd.hpp"

namespace cfd
{

template<typename T>
class finite_element
{
public:
  typedef T cell_type;
  static const unsigned int dimension = cell_type::dimension;
  typedef vertex<dimension> vertex_type;

  virtual evaluated_basis evaluate_basis(const std::size_t cid, const unsigned int i, const vertex_type& v) const = 0;
  virtual unsigned space_dimension() const = 0; // Number of basis functions
  virtual std::vector< std::pair<unsigned, unsigned> > getCommonDegreesOfFreedom(const cell_id cid, const cell_id cid2) const = 0;
  virtual std::vector<unsigned> getBoundaryDegreesOfFreedom(const cell_id cid, const std::vector< std::pair<vertex_id, vertex_id> >& boundary) const = 0;
  virtual ~finite_element() {}
};

}

#endif
