#include <boost/foreach.hpp>
#include <simple_cfd/symbolic/abstract_basic.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>

namespace cfd
{

namespace symbolic
{

namespace detail
{

Expr AbstractBasicHelper::integrate(const Expr& e, const Expr::region_t& region, const unsigned flags)
{
  Expr integrated = e;
  BOOST_FOREACH(const Expr::region_t::value_type& interval, region)
  {
    const Symbol& variable = interval.first;

    integrated = integrated.integrate(interval.first);
    Expr::subst_map lower, upper;
    lower[variable] = interval.second.first;
    upper[variable] = interval.second.second;

    integrated = (integrated.subs(upper) - integrated.subs(lower)).simplify();
  }

  return integrated;
}

}

}

}
