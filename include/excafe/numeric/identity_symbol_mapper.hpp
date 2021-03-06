#ifndef EXCAFE_NUMERIC_IDENTITY_SYMBOL_MAPPER_HPP
#define EXCAFE_NUMERIC_IDENTITY_SYMBOL_MAPPER_HPP

namespace excafe
{

namespace detail
{

template<typename T>
class IdentitySymbolMapper : public SymbolMapper<T, T>
{
public:
  typedef T value_type;

  value_type getSymbol(const value_type& k)
  {
    return k;
  }

  value_type getKey(const value_type& v) const
  {
    return v;
  }
};

}

}

#endif
