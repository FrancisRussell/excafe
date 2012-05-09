#ifndef EXCAFE_CAPTURE_FORMS_BILINEAR_FORM_INTEGRAL_SUM_HPP
#define EXCAFE_CAPTURE_FORMS_BILINEAR_FORM_INTEGRAL_SUM_HPP

#include <cassert>
#include <vector>
#include <utility>
#include <boost/shared_ptr.hpp>
#include <boost/variant/static_visitor.hpp>
#include <boost/variant/apply_visitor.hpp>
#include "field_expr.hpp"
#include "bilinear_form.hpp"
#include "bilinear_form_integral.hpp"

namespace excafe
{

namespace forms
{

class BilinearFormIntegralSum
{
private: 
  std::vector<BilinearForm> dxForms;
  std::vector<BilinearForm> dsForms;
  std::vector<BilinearForm> dSForms;

  class IntegralAdderVisitor : public boost::static_visitor<void>
  {
  private:
    BilinearFormIntegralSum& sum;
    const BilinearFormIntegral& integral;

  public:
    IntegralAdderVisitor(BilinearFormIntegralSum& _sum, const BilinearFormIntegral& i) : sum(_sum), integral(i)
    {
    }

    void operator()(const BilinearFormIntegral::cell_integral_tag&) const
    {
      sum.dxForms.push_back(integral.getForm());
    }

    void operator()(const BilinearFormIntegral::exterior_facet_integral_tag&) const
    {
      sum.dsForms.push_back(integral.getForm());
    }
  };

public:
  typedef std::vector<BilinearForm>::iterator iterator;
  typedef std::vector<BilinearForm>::const_iterator const_iterator;

  BilinearFormIntegralSum(const BilinearFormIntegral& i)
  {
    const BilinearFormIntegral::region_t region = i.getRegion();
    boost::apply_visitor(IntegralAdderVisitor(*this, i), region);
  }

  void append(const BilinearFormIntegralSum& b)
  {
    dxForms.insert(dxForms.end(), b.dxForms.begin(), b.dxForms.end());
    dsForms.insert(dsForms.end(), b.dsForms.begin(), b.dsForms.end());
    dSForms.insert(dsForms.end(), b.dSForms.begin(), b.dSForms.end());
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

  iterator begin_dS()
  {
    return dSForms.begin();
  }

  iterator end_dS()
  {
    return dSForms.end();
  }

  const_iterator begin_dS() const
  {
    return dSForms.begin();
  }

  const_iterator end_dS() const
  {
    return dSForms.end();
  }
};

}

}

#endif
