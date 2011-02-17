#ifndef SIMPLE_CFD_SYMBOLIC_ABSTRACT_BASIC_HPP
#define SIMPLE_CFD_SYMBOLIC_ABSTRACT_BASIC_HPP

#include <boost/functional/hash.hpp>
#include <boost/shared_ptr.hpp>
#include <simple_cfd/util/type_info.hpp>
#include "basic.hpp"
#include "visitor.hpp"

namespace cfd
{

namespace symbolic
{

template<typename T>
class AbstractBasic : public Basic
{
private:
  mutable bool isHashed;
  mutable std::size_t hash;
  bool heapAllocated;

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

  bool operator==(const AbstractBasic& b) const
  {
    return true;
  }

  AbstractBasic& operator=(const AbstractBasic& a)
  {
    isHashed = a.isHashed;
    hash = a.hash;
    return *this;
  }

  void invalidateHash()
  {
    isHashed = false;
  }

public:
  class Visitor
  {
  public:
    virtual void visit(const child_type& c) = 0;
  };

  bool operator==(const Basic& b) const
  {
    if (this == &b)
    {
      return true;
    }
    else if (hashValue() != b.hashValue())
    {
      return false;
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

  AbstractBasic() : isHashed(false), heapAllocated(false)
  {
  }

  AbstractBasic(const AbstractBasic& a) : 
    isHashed(a.isHashed), hash(a.hash), heapAllocated(false)
  {
  }

  void markHeapAllocated()
  {
    heapAllocated = true;
  }

  Expr clone() const
  {
    if (heapAllocated)
    {
      Expr::ref_t ref = this->shared_from_this();
      return Expr(ref);
    }
    else
    {
      child_type* const cloned = new child_type(asChild(*this));
      return Expr(cloned);
    }
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
    if (!isHashed)
    {
      hash = 0;
      boost::hash_combine(hash, getTypeHash());
      boost::hash_combine(hash, asChild(*this).untypedHash());
      isHashed=true;
    }
    return hash;
  }

  Expr expand() const
  {
    return clone();
  }

  void accept(symbolic::Visitor& v) const
  {
    if (Visitor* s = dynamic_cast<Visitor*>(&v))
      s->visit(asChild(*this));
    else
      v.visit(*this);
  }
};

}

}

#endif
