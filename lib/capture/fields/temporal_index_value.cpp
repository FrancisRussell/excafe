#include <cassert>
#include <simple_cfd/capture/fields/temporal_index_value.hpp>

namespace cfd
{

namespace detail
{

TemporalIndexValue::TemporalIndexValue() : terminationSet(false)
{
}

void TemporalIndexValue::setTermination(const ScalarExpr::expr_ptr c)
{
  assert(!terminationSet);
  terminationSet = true;
  termination = c;
}

ScalarExpr& TemporalIndexValue::getTermination() const
{
  return *termination;
}

void TemporalIndexValue::registerIndexable(IndexableValue<discrete_scalar_tag>& i)
{
  indexedScalars.insert(&i);
}

void TemporalIndexValue::registerIndexable(IndexableValue<discrete_field_tag>& i)
{
  indexedFields.insert(&i);
}

void TemporalIndexValue::registerIndexable(IndexableValue<discrete_operator_tag>& i)
{
  indexedOperators.insert(&i);
}

void TemporalIndexValue::unregisterIndexable(IndexableValue<discrete_scalar_tag>& i)
{
  indexedScalars.erase(&i);
}

void TemporalIndexValue::unregisterIndexable(IndexableValue<discrete_field_tag>& i)
{
  indexedFields.erase(&i);
}

void TemporalIndexValue::unregisterIndexable(IndexableValue<discrete_operator_tag>& i)
{
  indexedOperators.erase(&i);
}

std::set<IndexableValue<discrete_scalar_tag>*> TemporalIndexValue::getIndexableScalars() const
{
  return indexedScalars;
}

std::set<IndexableValue<discrete_field_tag>*> TemporalIndexValue::getIndexableFields() const
{
  return indexedFields;
}

std::set<IndexableValue<discrete_operator_tag>*> TemporalIndexValue::getIndexableOperators() const
{
  return indexedOperators;
}

}

}
