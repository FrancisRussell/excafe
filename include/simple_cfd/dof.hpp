#ifndef SIMPLE_CFD_DOF_HPP
#define SIMPLE_CFD_DOF_HPP

#include <ostream>
#include "simple_cfd_fwd.hpp"

namespace cfd
{

template<std::size_t D>
class Dof
{
private:
  static const std::size_t dimension = D;
  const FiniteElement<dimension>* element;
  std::size_t cid;
  std::size_t dof;
  
public:
  Dof(const FiniteElement<dimension>* const _element, const std::size_t _cid, const std::size_t _dof) : 
    element(_element), cid(_cid), dof(_dof)
  {
  }

  const FiniteElement<dimension>* getElement() const
  {
    return element;
  }

  std::size_t getCell() const
  {
    return cid;
  }

  std::size_t getIndex() const
  {
    return dof;
  }

  bool operator==(const Dof& d) const
  {
    return element == d.element && cid == d.cid && dof == d.dof;
  }

  bool operator<(const Dof& d) const
  {
    if (element < d.element) return true;
    if (element == d.element && cid < d.cid) return true;
    if (element == d.element && cid == d.cid && dof < d.dof) return true;
    return false;
  }

  void write(std::ostream& o) const
  {
    o << "dof(element=" << *element << ", cid=" << cid << ", index=" << dof << ")";
  }
};

}

namespace std
{

template<std::size_t D>
std::ostream& operator<<(std::ostream& o, const cfd::Dof<D>& dof)
{
  dof.write(o);
  return o;
}

}

#endif
