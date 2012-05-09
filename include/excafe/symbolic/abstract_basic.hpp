#ifndef EXCAFE_SYMBOLIC_ABSTRACT_BASIC_HPP
#define EXCAFE_SYMBOLIC_ABSTRACT_BASIC_HPP

#include <utility>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <excafe/util/type_info.hpp>
#include <excafe/util/hash.hpp>
#include "symbolic_fwd.hpp"
#include "basic.hpp"
#include "visitor.hpp"
#include "expr.hpp"
#include "make_expr_from.hpp"

namespace excafe
{

namespace symbolic
{

namespace detail
{

class AbstractBasicHelper
{
public:
  static Expr integrate(const Expr& e, const Expr::region_t& region, const unsigned flags);
};

}

template<typename T>
class AbstractBasic : public Basic
{
private:
  static const unsigned HEAP_ALLOCATED = 0x01;
  static const unsigned HASH_CACHED    = 0x02;

  mutable unsigned flags;
  mutable std::size_t hash;

protected:
  typedef T child_type;

  static util::TypeInfo getType(const Basic& b)
  {
    return util::TypeInfo(typeid(b));
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
    hash = a.hash;
    flags = (a.flags & ~HEAP_ALLOCATED) | (flags & HEAP_ALLOCATED);
    return *this;
  }

  void invalidateHash()
  {
    flags &= ~HASH_CACHED;
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
    else if (typeHash() != b.typeHash())
    {
      return false;
    }
    else if (hashValue() != b.hashValue())
    {
      return false;
    }
    else if (is_exactly_a<child_type>(b))
    {
      return asChild(*this)==asChild(b);
    }
    else
    {
      return false;
    }
  }

  AbstractBasic() : flags(0)
  {
  }

  AbstractBasic(const AbstractBasic& a) : 
    flags(a.flags & ~HEAP_ALLOCATED), hash(a.hash)
  {
  }

  std::size_t typeHash() const
  {
    static std::size_t typeHash = util::TypeInfo(typeid(child_type)).hashValue();
    return typeHash;
  }

  void markHeapAllocated()
  {
    flags |= HEAP_ALLOCATED;
  }

  Expr clone() const
  {
    if (flags & HEAP_ALLOCATED)
    {
      Expr::ref_t ref = this->shared_from_this();
      return Expr(ref);
    }
    else
    {
      return make_expr_from(asChild(*this));
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
    if (~flags & HASH_CACHED)
    {
      hash = 0x02c3866e;
      excafe::util::hash_accum(hash, typeHash());
      excafe::util::hash_accum(hash, asChild(*this).untypedHash());
      flags |= HASH_CACHED;
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

  Expr extractMultiplier(Rational& r) const
  {
    return this->simplify();
  }

  Expr integrate(const Expr::region_t& region, const unsigned flags) const
  {
    return detail::AbstractBasicHelper::integrate(clone(), region, flags);
  }
};

}

}

#endif
