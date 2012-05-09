#ifndef EXCAFE_CAPTURE_INDICES_PROPAGATION_RULES_HPP
#define EXCAFE_CAPTURE_INDICES_PROPAGATION_RULES_HPP

#include <set>
#include <map>
#include <memory>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include "propagation_rule.hpp"
#include "index_propagation_all.hpp"
#include "index_propagation_except.hpp"
#include "propagation_rule_visitor.hpp"
#include <excafe/capture/fields/temporal_index_set.hpp>
#include <excafe/capture/fields/discrete_expr.hpp>

namespace excafe
{

namespace detail
{

class RulePropagationHelper : public PropagationRuleVisitor
{
private:
  std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices;

public:
  RulePropagationHelper(std::map<DiscreteExpr*, TemporalIndexSet>& _exprIndices) : exprIndices(_exprIndices)
  {
  }

  void visit(IndexPropagationAll& i);
  void visit(IndexPropagationExcept& i);
};

class PropagationRules
{
private:
  typedef std::set< boost::shared_ptr<PropagationRule> > rule_set_t;
  rule_set_t rules;

public:
  typedef boost::indirect_iterator<rule_set_t::iterator> iterator;
  typedef boost::indirect_iterator<rule_set_t::const_iterator> const_iterator;

  void insert(std::auto_ptr<PropagationRule> rule)
  {
    rules.insert(boost::shared_ptr<PropagationRule>(rule.release()));
  }

  void propagateIndices(std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices) const;
  PropagationRules& operator+=(const PropagationRules& p);

  std::size_t size() const
  {
    return rules.size();
  }

  iterator begin()
  {
    return rules.begin();
  }

  iterator end()
  {
    return rules.end();
  }

  const_iterator begin() const
  {
    return rules.begin();
  }

  const_iterator end() const
  {
    return rules.end();
  }
};

}

}

#endif
