#ifndef SIMPLE_CFD_SYMBOLIC_DERIVATIVE_CACHE_HPP
#define SIMPLE_CFD_SYMBOLIC_DERIVATIVE_CACHE_HPP

#include "symbolic_fwd.hpp"
#include <utility>
#include <boost/scoped_ptr.hpp>
#include <boost/unordered_map.hpp>
#include <boost/utility.hpp>

namespace cfd
{

namespace symbolic
{

class DerivativeCache : public boost::noncopyable
{
private:
  typedef boost::unordered_map<std::pair<Expr, Symbol>, Expr> cache_t;
  boost::scoped_ptr<cache_t> cache;

public:
  DerivativeCache();
  Expr derivative(const Expr& e, const Symbol& s);
};

}

}

#endif
