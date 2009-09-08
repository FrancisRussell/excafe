#ifndef SIMPLE_CFD_CAPTURE_FIELDS_INDEXABLE_VALUE_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_INDEXABLE_VALUE_HPP

#include "temporal_index_expr.hpp"
#include "temporal_index.hpp"
#include "discrete_traits.hpp"
#include <map>
#include <cassert>
#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>

namespace cfd
{

namespace detail
{

template<typename discrete_object_tag>
class IndexableValue
{
public:
  typedef boost::shared_ptr<IndexableValue> value_ptr;
  typedef typename DiscreteTraits<discrete_object_tag>::expr_ptr expr_ptr;

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
      assert(offsetValue < 0 && "Can only assign variables for absolute index values below 0");
      assert(parent.initialValues.find(offsetValue) == parent.initialValues.end() && "Can only assign absolute value once");
      parent.initialValues.insert(std::make_pair(offsetValue, rhs)); 
    }

    void operator()(const TemporalIndexOffset::relative_tag&) const
    {
      assert(offsetValue == 0 && "Can only assign current value of an iteration");
      parent.assignedValue = rhs;
    }

    void operator()(const TemporalIndexOffset::final_tag&) const
    {
      assert(false && "Cannot assign to indexed value using relative-to-final index");
    }
  };

public:
  IndexableValue(const TemporalIndex _indexVariable) : indexVariable(_indexVariable.getIndex()),
    assignedValue(new typename DiscreteTraits<discrete_object_tag>::undefined_t())
  {
  }

  void handleAssignment(const TemporalIndexExpr& indexExpr, const expr_ptr& rhs)
  {
    assert(indexExpr.getIndex() == indexVariable && "Incorrect index used in expression lhs");
    
    TemporalIndexOffset offset = indexExpr.getOffset();
    TemporalIndexOffset::offset_t offsetType = offset.getType();

    OffsetTypeVisitor visitor(*this, offset.getValue(), rhs);
    boost::apply_visitor(visitor, offsetType);
  }

  expr_ptr getIterationAssignment() const
  {
    return assignedValue;
  }

  TemporalIndexValue::index_ptr getIndexVariable() const
  {
    return indexVariable;
  }
};

}

}

#endif
