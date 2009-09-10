#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SET_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_DISCRETE_EXPR_SET_HPP

#include <set>
#include <simple_cfd/capture/fields/discrete_traits.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class DiscreteExprSet
{
private:
  std::set<typename DiscreteTraits<discrete_object_tag>::expr_t*> expressions;
  std::set<typename DiscreteTraits<discrete_object_tag>::indexable_t*> indexedParents;

public:
  void addExpr(typename DiscreteTraits<discrete_object_tag>::expr_t& e)
  {
    expressions.insert(&e);
  }

  void addIndexed(typename DiscreteTraits<discrete_object_tag>::indexable_t& i)
  {
    indexedParents.insert(&i);
  }
};

}

}

#endif
