#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_SUM_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_SUM_HPP

#include <vector>
#include <utility>
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace forms
{

class BilinearFormSum
{
private:
  std::vector<std::pair<boost::shared_ptr<Field>, boost::shared_ptr<Field> > > forms;

public:
  BilinearFormSum(const BilinearForm& f)
  {
    forms.push_back(std::make_pair(f.getTrialField(), f.getTestField()));
  }

  void append(const BilinearFormSum& b)
  {
    forms.insert(forms.end(), b.forms.begin(), b.forms.end());
  }
};

}

}

#endif
