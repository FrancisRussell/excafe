#ifndef EXCAFE_UTIL_SINGLETON_HPP
#define EXCAFE_UTIL_SINGLETON_HPP

// Design inspired by boost/pool/detail/singleton.hpp

namespace excafe
{

namespace util
{

template<typename T>
class Singleton
{
public:
  typedef T value_type;

private:
  struct InstanceCreator
  {
    InstanceCreator() { value_type::getInstance(); }
    void forceInstantiation() {}
  };

  static InstanceCreator instanceCreator;

public:
  static value_type& getInstance()
  {
    instanceCreator.forceInstantiation();
    static value_type instance;
    return instance;
  }
};

}

}

#endif
