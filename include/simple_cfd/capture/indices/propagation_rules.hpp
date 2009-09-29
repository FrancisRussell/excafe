#ifndef SIMPLE_CFD_CAPTURE_INDICES_PROPOGATION_RULES_HPP
#define SIMPLE_CFD_CAPTURE_INDICES_PROPOGATION_RULES_HPP

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

namespace cfd
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

  void visit(IndexPropagationAll& i)
  {
    const std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator fromIter(exprIndices.find(&i.getFrom()));
    const std::map<DiscreteExpr*, TemporalIndexSet>::iterator toIter(exprIndices.find(&i.getTo()));

    assert(fromIter != exprIndices.end());
    assert(toIter != exprIndices.end());

    toIter->second += fromIter->second;
  }

  void visit(IndexPropagationExcept& i)
  {
    const std::map<DiscreteExpr*, TemporalIndexSet>::const_iterator fromIter(exprIndices.find(&i.getFrom()));
    const std::map<DiscreteExpr*, TemporalIndexSet>::iterator toIter(exprIndices.find(&i.getTo()));

    assert(fromIter != exprIndices.end());
    assert(toIter != exprIndices.end());

    TemporalIndexSet indices(fromIter->second);
    indices -= &i.getExcludedIndex();
    toIter->second += indices;
  }
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

  void propagateIndices(std::map<DiscreteExpr*, TemporalIndexSet>& exprIndices) const
  {
    bool converged = false;
    RulePropagationHelper propagationHelper(exprIndices);

    while(!converged)
    {
      const std::map<DiscreteExpr*, TemporalIndexSet> exprIndicesOld(exprIndices);

      for(const_iterator ruleIter(begin()); ruleIter!=end(); ++ruleIter)
      {
        ruleIter->accept(propagationHelper);
      }

      converged = (exprIndices == exprIndicesOld);
    }
  }

  PropagationRules& operator+=(const PropagationRules& p)
  {
    rules.insert(p.rules.begin(), p.rules.end());
    return *this;
  }

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
