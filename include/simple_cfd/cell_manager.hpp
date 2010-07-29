#ifndef SIMPLE_CFD_CELL_MANAGER_HPP
#define SIMPLE_CFD_CELL_MANAGER_HPP

#include <cstddef>
#include <simple_cfd/simple_cfd_fwd.hpp>
#include <simple_cfd/util/singleton.hpp>

namespace cfd
{

class CellManager
{
public:
  template<std::size_t D>
  struct ref { typedef const GeneralCell<D>* general; };

  typedef const MeshCell* mesh_cell_ref;

  template<typename cell_type>
  static typename ref<cell_type::dimension>::general getInstance()
  {
    return &(util::Singleton<cell_type>::getInstance());
  }
};

}

#endif
