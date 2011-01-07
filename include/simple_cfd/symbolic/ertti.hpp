#ifndef SIMPLE_CFD_ERTTI_HPP
#define SIMPLE_CFD_ERTTI_HPP

#include "type_manager.hpp"

namespace cfd
{

namespace symbolic
{

template<typename T, typename P>
class ERTTI
{
private:
  typedef T child_type;

  struct TypeIDInitialiser
  {
    TypeIDInitialiser()
    {
      ERTTI<T,P>::typeID();
    }

    void doNothing() {}
  };

  static TypeIDInitialiser initialiser;


public:
  static int typeID()
  {
    initialiser.doNothing();
    static int id = TypeManager::getNewTypeID();
    //static int parentId = P::typeID();
    return id;
  };

  virtual int getTypeID() const
  {
    return typeID();
  }
};

}

}

#endif
