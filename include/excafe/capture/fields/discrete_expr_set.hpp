#ifndef EXCAFE_CAPTURE_EVALUATION_DISCRETE_EXPR_SET_HPP
#define EXCAFE_CAPTURE_EVALUATION_DISCRETE_EXPR_SET_HPP

#include <set>
#include <excafe/capture/fields/discrete_traits.hpp>
#include <boost/iterator/indirect_iterator.hpp>

namespace excafe
{

namespace detail
{

template<typename discrete_object_tag>
class DiscreteExprSet
{
private:
  typedef std::set<typename DiscreteTraits<discrete_object_tag>::expr_t*> expr_set_t;
  typedef std::set<typename DiscreteTraits<discrete_object_tag>::indexable_t*> indexable_set_t;

  expr_set_t expressions;
  indexable_set_t indexedParents;

public:
  typedef boost::indirect_iterator<typename expr_set_t::iterator> expr_iter;
  typedef boost::indirect_iterator<typename indexable_set_t::iterator> indexable_iter;

  bool insert(typename DiscreteTraits<discrete_object_tag>::expr_t& e)
  {
    return expressions.insert(&e).second;
  }

  bool insert(typename DiscreteTraits<discrete_object_tag>::indexable_t& i)
  {
    return indexedParents.insert(&i).second;
  }

  expr_iter begin_expr() const
  {
    return expr_iter(expressions.begin());
  }

  expr_iter end_expr() const
  {
    return expr_iter(expressions.end());
  }

  indexable_iter begin_indexable() const
  {
    return indexable_iter(indexedParents.begin());
  }

  indexable_iter end_indexable() const
  {
    return indexable_iter(indexedParents.end());
  }
};

}

}

#endif
