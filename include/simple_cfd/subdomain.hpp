#ifndef SIMPLE_CFD_SUBDOMAIN_HPP
#define SIMPLE_CFD_SUBDOMAIN_HPP

#include "simple_cfd_fwd.hpp"

namespace cfd
{

template<unsigned int D>
class SubDomain
{
public:
  static const unsigned int dimension = D;

  virtual bool inside(const vertex<dimension>& v) = 0;
  virtual ~SubDomain() {}
};

}

#endif
