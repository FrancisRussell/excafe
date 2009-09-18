#ifndef SIMPLE_CFD_CAPTURE_FIELDS_INDEXABLE_VALUE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_INDEXABLE_VALUE_HPP

#include <map>
#include <utility>
#include <iterator>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "temporal_index_set.hpp"
#include "temporal_index_expr.hpp"
#include "temporal_index.hpp"
#include "discrete_traits.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class IndexableValueInitialisationIterator : public boost::iterator_facade<
  IndexableValueInitialisationIterator<discrete_object_tag>, 
  typename DiscreteTraits<discrete_object_tag>::expr_t,
  boost::bidirectional_traversal_tag>
{
private:
  typedef typename DiscreteTraits<discrete_object_tag>::expr_t expr_t;
  typedef typename DiscreteTraits<discrete_object_tag>::expr_ptr expr_ptr;
  typedef typename std::map<signed, expr_ptr>::const_iterator iterator_t;
  iterator_t pos;

public:
  IndexableValueInitialisationIterator(const iterator_t& _pos) : pos(_pos) 
  {
  }

  signed getOffset() const
  {
    return pos->first;
  }

  expr_t& dereference() const
  {
    return *pos->second;
  }

  void increment()
  {
    ++pos;
  }

  void decrement()
  {
    --pos;
  }

  bool equal(const IndexableValueInitialisationIterator& i) const
  {
    return pos == i.pos;
  }
};

template<typename discrete_object_tag>
class IndexableValue
{
public:
  typedef boost::shared_ptr<IndexableValue> value_ptr;
  typedef typename DiscreteTraits<discrete_object_tag>::expr_ptr expr_ptr;
  typedef IndexableValueInitialisationIterator<discrete_object_tag> init_iterator;

private:
  const TemporalIndexValue::index_ptr indexVariable;
  std::map<signed, expr_ptr> initialValues;
  expr_ptr assignedValue;

  class OffsetTypeVisitor : public boost::static_visitor<void>
  {
  private:
    IndexableValue& parent;
    const signed offsetValue;
    const expr_ptr rhs;

  public:
    OffsetTypeVisitor(IndexableValue& _parent, const signed _offsetValue, const expr_ptr _rhs) : 
      parent(_parent), offsetValue(_offsetValue), rhs(_rhs)
    {
    }

    void operator()(const TemporalIndexOffset::absolute_tag&) const
    {
      if (offsetValue >= 0)
        CFD_EXCEPTION("Can only assign variables for absolute index values below 0");

      if (parent.initialValues.find(offsetValue) != parent.initialValues.end())
        CFD_EXCEPTION("Can only assign absolute value once");

      parent.initialValues.insert(std::make_pair(offsetValue, rhs)); 
    }

    void operator()(const TemporalIndexOffset::relative_tag&) const
    {
      if (offsetValue != 0)
        CFD_EXCEPTION("Can only assign current value of an iteration");

      parent.assignedValue = rhs;
    }

    void operator()(const TemporalIndexOffset::final_tag&) const
    {
      CFD_EXCEPTION("Cannot assign to indexed value using relative-to-final index");
    }
  };

public:
  IndexableValue(const TemporalIndex _indexVariable) : indexVariable(_indexVariable.getIndex()),
    assignedValue(new typename DiscreteTraits<discrete_object_tag>::undefined_t())
  {
  }

  void handleAssignment(const TemporalIndexExpr& indexExpr, const expr_ptr& rhs)
  {
    if(indexExpr.getIndex() != indexVariable)
      CFD_EXCEPTION("Incorrect index used in expression lhs");
    
    TemporalIndexOffset offset = indexExpr.getOffset();
    TemporalIndexOffset::offset_t offsetType = offset.getType();

    OffsetTypeVisitor visitor(*this, offset.getValue(), rhs);
    boost::apply_visitor(visitor, offsetType);
  }

  init_iterator begin_inits() const
  {
    return init_iterator(initialValues.begin());
  }

  init_iterator end_inits() const
  {
    return init_iterator(initialValues.end());
  }

  expr_ptr getIterationAssignment() const
  {
    return assignedValue;
  }

  TemporalIndexValue::index_ptr getIndexVariable() const
  {
    return indexVariable;
  }

  TemporalIndexSet getTemporalIndices() const
  {
    TemporalIndexSet indices;

    for(init_iterator initIter(begin_inits()); initIter!=end_inits(); ++initIter)
    {
      indices += initIter->getTemporalIndices();
    }

    indices += assignedValue->getTemporalIndices();
    indices += &(*indexVariable);
    
    assert(indices.contains(indexVariable->getTermination().getTemporalIndices()));

    return indices;
  }
};

}

}

#endif
