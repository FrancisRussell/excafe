#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_SUM_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_SUM_HPP

#include <vector>
#include <utility>
#include <boost/shared_ptr.hpp>
#include "field.hpp"
#include "bilinear_form.hpp"

namespace cfd
{

namespace forms
{

class BilinearFormSum
{
private: 
  std::vector<BilinearForm> forms;

public:
  typedef std::vector<BilinearForm>::iterator iterator;
  typedef std::vector<BilinearForm>::const_iterator const_iterator;

  BilinearFormSum(const BilinearForm& f)
  {
    forms.push_back(f);
  }

  void append(const BilinearFormSum& b)
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
