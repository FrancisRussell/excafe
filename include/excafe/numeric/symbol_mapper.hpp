#ifndef EXCAFE_NUMERIC_SYMBOL_MAPPER_HPP
#define EXCAFE_NUMERIC_SYMBOL_MAPPER_HPP

namespace excafe
{

namespace detail
{

template<typename K, typename V>
class SymbolMapper
{
public:
  typedef K key_type;
  typedef V value_type;

  virtual value_type getSymbol(const key_type& k) = 0;
  virtual key_type   getKey(const value_type& v) const = 0;
};

}

}

#endif
