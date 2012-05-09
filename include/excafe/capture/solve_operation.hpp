#ifndef EXCAFE_CAPTURE_SOLVE_OPERATION_HPP
#define EXCAFE_CAPTURE_SOLVE_OPERATION_HPP

#include <map>
#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "capture_fwd.hpp"
#include "fields/fields_fwd.hpp"
#include "fields/named_field.hpp"
#include "fields/field.hpp"
#include "evaluation/evaluation_strategy.hpp"

namespace excafe
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
    detail::ExpressionValues<D> evaluatedValues = evaluationStrategy->execute<D>(scenario);

    for (std::map<NamedField, Field>::const_iterator newIter(newValues.begin()); newIter!=newValues.end(); ++newIter)
    {
      scenario.setNamedValue(newIter->first.getName(), evaluatedValues.getValue(*newIter->second.getExpr()));
    }
  }
};

}

#endif
