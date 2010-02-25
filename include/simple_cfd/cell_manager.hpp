#ifndef SIMPLE_CFD_CELL_MANAGER_HPP
#define SIMPLE_CFD_CELL_MANAGER_HPP

#include <cstddef>
#include <boost/scoped_ptr.hpp>
#include <simple_cfd/simple_cfd_fwd.hpp>

namespace cfd
{

class CellManager
{
private:
  template<typename T>
  class Singleton
  {
  private:
    typedef T cell_type;
    typedef boost::scoped_ptr< GeneralCell<cell_type::dimension> > holder_t;

  public:
    static GeneralCell<T::dimension>* getInstance()
    {
      static holder_t instance;
      if (instance == false)
      {
        holder_t created(new cell_type());
        instance.swap(created);
      }
      return instance.get(); 
    }
  };

public:
  template<std::size_t D>
  struct ref { typedef const GeneralCell<D>* general; };


  typedef const MeshCell* mesh_cell_ref;

  template<typename cell_type>
  static typename ref<cell_type::dimension>::general getInstance()
  {
    return Singleton<cell_type>::getInstance();
  }
};

}

#endif
