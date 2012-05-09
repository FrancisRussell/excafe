#ifndef EXCAFE_SUBDOMAIN_HPP
#define EXCAFE_SUBDOMAIN_HPP

#include "excafe_fwd.hpp"

namespace excafe
{

template<unsigned int D>
class SubDomain
{
public:
  static const unsigned int dimension = D;

  virtual bool inside(const vertex<dimension>& v) const = 0;
  virtual ~SubDomain() {}
};

}

#endif
