#ifndef EXCAFE_CAPTURE_FIELDS_TEMPORAL_INDEX_VALUE_HPP
#define EXCAFE_CAPTURE_FIELDS_TEMPORAL_INDEX_VALUE_HPP

#include <set>
#include "fields_fwd.hpp"
#include "scalar_expr.hpp"
#include "discrete_traits.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/utility.hpp>

namespace excafe
{

namespace detail
{

class TemporalIndexValue : public boost::noncopyable
{
private:
  bool terminationSet;
  ScalarExpr::expr_ptr termination;
  std::set<IndexableValue<discrete_scalar_tag>*>   indexedScalars;
  std::set<IndexableValue<discrete_field_tag>*>    indexedFields;
  std::set<IndexableValue<discrete_operator_tag>*> indexedOperators;

public:
  typedef boost::shared_ptr<TemporalIndexValue> index_ptr;

  TemporalIndexValue();
  void setTermination(const ScalarExpr::expr_ptr c);
  ScalarExpr& getTermination() const;

  void registerIndexable(IndexableValue<discrete_scalar_tag>& i);
  void registerIndexable(IndexableValue<discrete_field_tag>& i);
  void registerIndexable(IndexableValue<discrete_operator_tag>& i);

  void unregisterIndexable(IndexableValue<discrete_scalar_tag>& i);
  void unregisterIndexable(IndexableValue<discrete_field_tag>& i);
  void unregisterIndexable(IndexableValue<discrete_operator_tag>& i);

  std::set<IndexableValue<discrete_scalar_tag>*> getIndexableScalars() const;
  std::set<IndexableValue<discrete_field_tag>*> getIndexableFields() const;
  std::set<IndexableValue<discrete_operator_tag>*> getIndexableOperators() const;
};

}

}

#endif
