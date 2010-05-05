#ifndef SIMPLE_CFD_NUMERIC_GINAC_MAPPER_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_MAPPER_HPP

#include <map>
#include <sstream>
#include <ginac/symbol.h>

namespace cfd
{

namespace detail
{

template<typename T>
class GinacMapper
{
private:
  typedef T key_type;
  typedef GiNaC::symbol value_type;

  std::map<key_type, value_type> mappings;

  struct object_creator
  {
    object_creator() { GinacMapper<T>::instance(); }
    inline void do_nothing() const {}
  };
  static object_creator create_object;

  GinacMapper()
  {
  }

public:
  static GinacMapper& instance()
  {
    static GinacMapper mapper;
    create_object.do_nothing();
    return mapper;
  }

  value_type getSymbol(const key_type& k)
  {
    const typename std::map<key_type, value_type>::iterator iter = mappings.find(k);

    if (iter != mappings.end())
    {
      return iter->second;
    }
    else
    {
      std::ostringstream name;
      name << k;
      const value_type symbol(name.str());
      mappings.insert(std::make_pair(k, symbol));
      return symbol;
    }
  }
};

template<typename T>
typename GinacMapper<T>::object_creator GinacMapper<T>::create_object;

}

}

#endif
