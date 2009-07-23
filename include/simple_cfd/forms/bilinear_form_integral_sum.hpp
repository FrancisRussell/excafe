#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_INTEGRAL_SUM_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_INTEGRAL_SUM_HPP

#include <cassert>
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
  std::vector<BilinearForm> dxForms;
  std::vector<BilinearForm> dsForms;

public:
  typedef std::vector<BilinearForm>::iterator iterator;
  typedef std::vector<BilinearForm>::const_iterator const_iterator;

  BilinearFormIntegralSum(const BilinearFormIntegral& i)
  {
    if (i.isCellIntegral())
    {
      dxForms.push_back(i.getForm());
    }
    else if (i.isExteriorFacetIntegral())
    {
      dsForms.push_back(i.getForm());
    }
    else
    {
      assert(false);
    }
  }

  void append(const BilinearFormIntegralSum& b)
  {
    dxForms.insert(dxForms.end(), b.dxForms.begin(), b.dxForms.end());
    dsForms.insert(dsForms.end(), b.dsForms.begin(), b.dsForms.end());
  }

  iterator begin_dx()
  {
    return dxForms.begin();
  }

  iterator end_dx()
  {
    return dxForms.end();
  }

  const_iterator begin_dx() const
  {
    return dxForms.begin();
  }

  const_iterator end_dx() const
  {
    return dxForms.end();
  }

  iterator begin_ds()
  {
    return dsForms.begin();
  }

  iterator end_ds()
  {
    return dsForms.end();
  }

  const_iterator begin_ds() const
  {
    return dsForms.begin();
  }

  const_iterator end_ds() const
  {
    return dsForms.end();
  }
};

}

}

#endif
