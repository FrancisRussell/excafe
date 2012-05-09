#ifndef EXCAFE_NUMERIC_EXCAFE_MAPPER_HPP
#define EXCAFE_NUMERIC_EXCAFE_MAPPER_HPP

#include <sstream>
#include <ostream>
#include <map>
#include <utility>
#include <boost/noncopyable.hpp>
#include <excafe/exception.hpp>
#include <excafe/util/singleton.hpp>
#include <excafe/symbolic/symbol.hpp>
#include "symbol_mapper.hpp"

namespace excafe
{

namespace detail
{

template<typename T>
class ExcafeMapper : boost::noncopyable, public SymbolMapper<T, symbolic::Symbol>
{
private:
  friend class util::Singleton< ExcafeMapper<T> >;
  typedef T key_type;
  typedef symbolic::Symbol value_type;

  std::map<key_type, value_type> mappings;
  std::map<value_type, key_type> reverse_mappings;

  ExcafeMapper()
  {
  }

  void checkConsistent() const
  {
    if (mappings.size() != reverse_mappings.size())
      CFD_EXCEPTION("Inconsistency detected in ExcafeMapper bi-directional mapping.");
  }

public:
  static ExcafeMapper& instance()
  {
    return util::Singleton<ExcafeMapper>::getInstance();
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
      CFD_EXCEPTION("Attempted to convert an unknown Excafe symbol back to a key.");
    }
  }
};

}

}

#endif
