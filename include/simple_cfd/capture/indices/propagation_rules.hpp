#ifndef SIMPLE_CFD_CAPTURE_INDICES_PROPOGATION_RULES_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_PROPOGATION_RULES_HPP

#include <memory>
#include <boost/ptr_container/ptr_set.hpp>
#include "propagation_rule.hpp"

namespace cfd
{

namespace detail
{

class PropagationRules
{
private:
  boost::ptr_set<PropagationRule> rules;

public:
  void insert(std::auto_ptr<PropagationRule> rule)
  {
    rules.insert(rule.release());
  }

  PropagationRules& operator+=(const PropagationRules& p)
  {
    rules.insert(p.rules.begin(), p.rules.end());
    return *this;
  }
};

}

}

#endif
