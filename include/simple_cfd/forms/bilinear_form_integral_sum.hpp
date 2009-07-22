#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_INTEGRAL_SUM_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_INTEGRAL_SUM_HPP

#include <vector>
#include <utility>
#include <boost/shared_ptr.hpp>
#include "field.hpp"
#include "bilinear_form.hpp"
#include "bilinear_form_integral.hpp"

namespace cfd
{

namespace forms
{

class BilinearFormIntegralSum
{
private: 
  std::vector<BilinearForm> forms;

public:
  typedef std::vector<BilinearForm>::iterator iterator;
  typedef std::vector<BilinearForm>::const_iterator const_iterator;

  BilinearFormIntegralSum(const BilinearFormIntegral& i)
  {
    forms.push_back(i.getForm());
  }

  void append(const BilinearFormIntegralSum& b)
  {
    forms.insert(forms.end(), b.forms.begin(), b.forms.end());
  }

  iterator begin()
  {
    return forms.begin();
  }

  iterator end()
  {
    return forms.end();
  }

  const_iterator begin() const
  {
    return forms.begin();
  }

  const_iterator end() const
  {
    return forms.end();
  }
};

}

}

#endif
