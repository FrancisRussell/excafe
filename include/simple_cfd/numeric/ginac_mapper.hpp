#ifndef SIMPLE_CFD_NUMERIC_GINAC_MAPPER_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_MAPPER_HPP

#include <sstream>
#include <map>
#include <utility>
#include <boost/noncopyable.hpp>
#include <simple_cfd/exception.hpp>
#include <ginac/symbol.h>

namespace cfd
{

namespace detail
{

template<typename T>
class GinacMapper : boost::noncopyable
{
private:
  typedef T key_type;
  typedef GiNaC::symbol value_type;

  std::map<key_type, value_type> mappings;
  
  // NOTE: we need to use a GiNaC::ex_is_less to get a sensible ordering
  std::map<value_type, key_type, GiNaC::ex_is_less> reverse_mappings;

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

  value_type getGiNaCSymbol(const key_type& k)
  {
    const typename std::map<key_type, value_type>::const_iterator iter = mappings.find(k);

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
      reverse_mappings.insert(std::make_pair(symbol, k));
      return symbol;
    }
  }

  key_type getOriginalSymbol(const value_type& v) const
  {
    const typename std::map<value_type, key_type>::const_iterator iter = reverse_mappings.find(v);

    if (iter != reverse_mappings.end())
    {
      return iter->second;
    }
    else
    {
      CFD_EXCEPTION("Attempted to convert an unknown GiNaC symbol back to a key");
    }
  }
};

template<typename T>
typename GinacMapper<T>::object_creator GinacMapper<T>::create_object;

}

}

#endif
