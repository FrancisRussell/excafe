#ifndef SIMPLE_CFD_SYMBOLIC_TYPE_MANAGER_HPP
#define SIMPLE_CFD_SYMBOLIC_TYPE_MANAGER_HPP

#include <map>
#include <utility>
#include <boost/utility.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <simple_cfd/util/type_info.hpp>
#include <simple_cfd/util/singleton.hpp>

namespace cfd
{

namespace symbolic
{

class TypeManager : public boost::noncopyable
{
private:
  friend class util::Singleton<TypeManager>;

  typedef std::map<util::TypeInfo, std::size_t> id_map_t;
  id_map_t idMap;

  TypeManager()
  {
  }

public:
  static TypeManager& getInstance()
  {
    return util::Singleton<TypeManager>::getInstance();
  }

  template<typename T>
  std::size_t typeID()
  {
    BOOST_STATIC_ASSERT((boost::is_base_of<Basic,T>::value));

    const util::TypeInfo typeInfo(typeid(T));
    const id_map_t::const_iterator iter = idMap.find(typeInfo);
    {
      if (iter == idMap.end())
      {
        const std::size_t newID = idMap.size();
        idMap.insert(std::make_pair(typeInfo, newID));
        return newID;
      }
      else
      {
        return iter->second;
      }
    }
  }
};

}

}
#endif
