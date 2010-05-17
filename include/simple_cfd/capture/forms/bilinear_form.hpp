#ifndef SIMPLE_CFD_CAPTURE_FORMS_BILINEAR_FORM_HPP
#define SIMPLE_CFD_CAPTURE_FORMS_BILINEAR_FORM_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "linear_form.hpp"
#include "field_expr.hpp"

namespace cfd
{

namespace forms
{

class BilinearForm
{
private:
  boost::shared_ptr< detail::FieldExpr > trialField;
  boost::shared_ptr< detail::FieldExpr > testField;

public:

  BilinearForm(const LinearForm& trial, const LinearForm& test) :
    trialField(trial.getField()), testField(test.getField())
  {
  }

  boost::shared_ptr<detail::FieldExpr> getTrialField() const
  {
    return trialField;
  }

  boost::shared_ptr<detail::FieldExpr> getTestField() const
  {
    return testField;
  }
};

}

}

#endif
