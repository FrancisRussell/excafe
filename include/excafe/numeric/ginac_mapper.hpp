#ifndef EXCAFE_NUMERIC_GINAC_MAPPER_HPP
#define EXCAFE_NUMERIC_GINAC_MAPPER_HPP

#include <sstream>
#include <map>
#include <utility>
#include <boost/noncopyable.hpp>
#include <excafe/exception.hpp>
#include <excafe/util/singleton.hpp>
#include "symbol_mapper.hpp"
#include <ginac/symbol.h>

namespace excafe
{

namespace detail
{

template<typename T>
class GinacMapper : boost::noncopyable, public SymbolMapper<T, GiNaC::symbol>
{
private:
  friend class util::Singleton< GinacMapper<T> >;
  typedef T key_type;
  typedef GiNaC::symbol value_type;

  std::map<key_type, value_type> mappings;
  
  // NOTE: we need to use a GiNaC::ex_is_less to get a consistent ordering
  std::map<value_type, key_type, GiNaC::ex_is_less> reverse_mappings;

  GinacMapper()
  {
  }

  void checkConsistent() const
  {
    if (mappings.size() != reverse_mappings.size())
      CFD_EXCEPTION("Inconsistency detected in GinacMapper bi-directional mapping.");
  }

public:
  static GinacMapper& instance()
  {
    return util::Singleton<GinacMapper>::getInstance();
  }

  value_type getSymbol(const key_type& k)
  {
    checkConsistent();

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

  key_type getKey(const value_type& v) const
  {
    checkConsistent();

    const typename std::map<value_type, key_type>::const_iterator iter = reverse_mappings.find(v);

    if (iter != reverse_mappings.end())
    {
      return iter->second;
    }
    else
    {
      CFD_EXCEPTION("Attempted to convert an unknown GiNaC symbol back to a key.");
    }
  }
};

}

}

#endif
