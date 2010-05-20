#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_MAP_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_EXPRESSION_VALUES_MAP_HPP

#include <cassert>
#include <cstddef>
#include <map>
#include <utility>
#include <simple_cfd/discrete_value_traits.hpp>
#include <simple_cfd/capture/fields/temporal_index_expr.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag, std::size_t D>
class ExpressionValuesMap
{
private:
  static const std::size_t dimension = D;
  typedef typename DiscreteTraits<discrete_object_tag>::expr_t expr_t;
  typedef typename DiscreteTraits<discrete_object_tag>::indexable_t indexable_t;
  typedef typename DiscreteValueTraits<discrete_object_tag, dimension>::value_t value_t;

  std::map<expr_t*, value_t> nonIndexed;
  std::map<std::pair<indexable_t*, int>, value_t> indexed;

public:
  void completeIteration()
  {
    nonIndexed.clear();
    
    std::map<std::pair<indexable_t*, int>, value_t> newIndexed;

    for(typename std::map<std::pair<indexable_t*, int>, value_t>::const_iterator indexedIter(indexed.begin());
      indexedIter != indexed.end(); ++indexedIter)
    {
      const std::pair<indexable_t*, int> oldKey(indexedIter->first);
      const std::pair<indexable_t*, int> newKey(std::make_pair(oldKey.first, oldKey.second-1));
      newIndexed.insert(std::make_pair(newKey, indexedIter->second));
    }

    std::swap(indexed, newIndexed);
  }

  void addMappings(const std::map<expr_t*, value_t> newValues)
  {
    nonIndexed.insert(newValues.begin(), newValues.end());
  }

  bool hasValue(expr_t& e) const
  {
    return nonIndexed.find(&e) != nonIndexed.end();
  }

  bool hasValue(indexable_t& e, const signed offset) const
  {
    return indexed.find(std::make_pair(&e, offset)) != indexed.end();
  }
  
  value_t& getValue(expr_t& e)
  {
    const typename std::map<expr_t*, value_t>::iterator exprIter = nonIndexed.find(&e); 
    assert(exprIter != nonIndexed.end());
    return exprIter->second;
  }

  value_t& getValue(indexable_t& e, const signed offset)
  {
    const typename std::map<std::pair<indexable_t*, int>, value_t>::iterator exprIter = 
      indexed.find(std::make_pair(&e, offset));

    assert(exprIter != indexed.end());
    return exprIter->second;
  }

  void setValue(expr_t& e, const value_t& v)
  {
    assert(!hasValue(e));
    nonIndexed.insert(std::make_pair(&e, v));
  }
  
  void setValue(indexable_t& e, const value_t& v, const signed offset)
  {
    assert(!hasValue(e, offset));
    const std::pair<indexable_t*, int> key(&e, offset);
    indexed.insert(std::make_pair(key, v));
  }

  std::map<expr_t*, value_t> getFinals() const
  {
    std::map<expr_t*, value_t> finals;

    for(typename std::map<std::pair<indexable_t*, int>, value_t>::const_iterator indexedIter(indexed.begin());
      indexedIter != indexed.end(); ++indexedIter)
    {
      indexable_t& indexable = *indexedIter->first.first;
      const signed offset = indexedIter->first.second;

      const TemporalIndexExpr indexExpr(TemporalIndexExpr::final(indexable.getIndexVariable(), offset));
      const typename DiscreteTraits<discrete_object_tag>::expr_ptr expr = indexable.getIndexedExpr(indexExpr);

      if (expr.use_count() > 0)
      {
        finals.insert(std::make_pair(&(*expr), indexedIter->second));
      }
    }

    return finals;
  }
};

}

}

#endif
