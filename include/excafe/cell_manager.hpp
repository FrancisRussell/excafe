#ifndef EXCAFE_CELL_MANAGER_HPP
#define EXCAFE_CELL_MANAGER_HPP

#include <cstddef>
#include <excafe/excafe_fwd.hpp>
#include <excafe/util/singleton.hpp>

namespace excafe
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
