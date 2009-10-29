#ifndef SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP
#define SIMPLE_CFD_CAPTURE_SOLVE_OPERATION_HPP

#include <map>
#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "capture_fwd.hpp"
#include "fields/fields_fwd.hpp"
#include "fields/named_field.hpp"
#include "fields/field.hpp"
#include "evaluation/evaluation_strategy.hpp"

namespace cfd
{

class SolveOperation
{
private:
  detail::DimensionlessScenario* scenario;
  bool finished;
  std::map<NamedField, Field> newValues;

  //NOTE: we want ownership here, but have to support assignment
  boost::shared_ptr<detail::EvaluationStrategy> evaluationStrategy;

  bool isExecutable() const;

public:
  SolveOperation();
  SolveOperation(detail::DimensionlessScenario& s);
  void setNewValue(NamedField& f, const Field& value);
  void finish();
  void execute();

  template<std::size_t D>
  void executeDimensionTemplated(Scenario<D>& scenario)
  {
    evaluationStrategy->execute<D>(scenario);
  }
};

}

#endif
