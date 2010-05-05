#ifndef SIMPLE_CFD_NUMERIC_GINAC_MAPPER_HPP
#define SIMPLE_CFD_NUMERIC_GINAC_MAPPER_HPP

#include <sstream>
#include <boost/bimap.hpp>
#include <simple_cfd/exception.hpp>
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

  boost::bimap<key_type, value_type> mappings;

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
    const typename boost::bimap<key_type, value_type>::left_iterator iter = mappings.left.find(k);

    if (iter != mappings.left.end())
    {
      return iter->second;
    }
    else
    {
      std::ostringstream name;
      name << k;
      const value_type symbol(name.str());
      mappings.left.insert(std::make_pair(k, symbol));
      return symbol;
    }
  }

  key_type getOriginalSymbol(const value_type& v) const
  {
    const typename boost::bimap<key_type, value_type>::right_const_iterator iter = mappings.right.find(v);

    if (iter != mappings.right.end())
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
