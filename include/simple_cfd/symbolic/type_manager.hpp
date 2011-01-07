#ifndef SIMPLE_CFD_SYMBOLIC_TYPE_MANAGER_HPP
#define SIMPLE_CFD_SYMBOLIC_TYPE_MANAGER_HPP

namespace cfd
{

namespace symbolic
{

class TypeManager
{
public:
  static int getNewTypeID()
  {
    // If this is a static member, we can't be certain that it's been 
    // initialised when this function is called
    static int nextID = 0;
    return nextID++;
  }
};

}

}

#endif
