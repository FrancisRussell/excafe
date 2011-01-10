#ifndef SIMPLE_CFD_SYMBOLIC_ABSTRACT_BASIC_HPP
#define SIMPLE_CFD_SYMBOLIC_ABSTRACT_BASIC_HPP

#include <boost/functional/hash.hpp>
#include <simple_cfd/util/type_info.hpp>
#include "basic.hpp"

namespace cfd
{

namespace symbolic
{

template<typename T>
class AbstractBasic : public Basic
{
protected:
  typedef T child_type;

  static util::TypeInfo getType(const Basic& b)
  {
    return util::TypeInfo(typeid(b));
  }

  static std::size_t getTypeHash()
  {
    static std::size_t typeHash = util::TypeInfo(typeid(child_type)).hashValue();
    return typeHash;
  }

  static const child_type& asChild(const Basic& b)
  {
    return static_cast<const child_type&>(b);
  }
  
  static child_type& asChild(Basic& b)
  {
    return static_cast<child_type&>(b);
  }

public:
  bool operator==(const Basic& b) const
  {
    if (this == &b)
    {
      return true;
    }
    else if (getType(*this) == getType(b))
    {
      return asChild(*this)==asChild(b);
    }
    else
    {
      return false;
    }
  }

  bool operator<(const Basic& b) const
  {
    if (this == &b)
    {
      return false;
    }
    else if (getType(*this) < getType(b))
    {
      return true;
    }
    else if (getType(*this) == getType(b))
    {
      return asChild(*this) < asChild(b);
    }
    else
    {
      return false;
    }
  }

  Expr clone() const
  {
    return Expr(new child_type(asChild(*this)));
  }

  Expr simplify() const
  {
    return clone();
  }

  operator Expr() const
  {
    return clone();
  }

  std::size_t hashValue() const
  {
    std::size_t hash = getTypeHash();
    boost::hash_combine(hash, asChild(*this).untypedHash());
    return hash;
  }
};

}

}

#endif
