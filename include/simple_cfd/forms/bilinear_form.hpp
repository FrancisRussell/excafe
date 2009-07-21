#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_HPP

#include <cstddef>
#include <boost/shared_ptr.hpp>
#include "linear_form.hpp"
#include "field.hpp"

namespace cfd
{

namespace forms
{

class BilinearForm
{
private:
  boost::shared_ptr< Field > trialField;
  boost::shared_ptr< Field > testField;

public:

  BilinearForm(const LinearForm& trial, const LinearForm& test) :
    trialField(trial.getField()), testField(test.getField())
  {
  }

  boost::shared_ptr<Field> getTrialField() const
  {
    return trialField;
  }

  boost::shared_ptr<Field> getTestField() const
  {
    return testField;
  }
};

}

}

#endif
