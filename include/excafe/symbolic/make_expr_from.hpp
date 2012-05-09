#ifndef EXCAFE_SYMBOLIC_MAKE_EXPR_FROM_HPP
#define EXCAFE_SYMBOLIC_MAKE_EXPR_FROM_HPP

#include "basic.hpp"
#include "expr.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace excafe
{

namespace symbolic
{

template<typename T>
Expr make_expr_from(const T& value)
{
  const boost::shared_ptr<T> ptr = boost::make_shared<T>(value);
  ptr->markHeapAllocated();
  const Expr::ref_t ref(ptr);
  return Expr(ref);
}

}

}

#endif
