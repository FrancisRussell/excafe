#ifndef SIMPLE_CFD_FORMS_FORMS_HPP
#define SIMPLE_CFD_FORMS_FORMS_HPP

#include <cstddef>
#include <cassert>
#include <boost/shared_ptr.hpp>

#include "field.hpp"
#include "linear_form.hpp"
#include "bilinear_form.hpp"
#include "bilinear_form_integral.hpp"
#include "bilinear_form_integral_sum.hpp"
#include "facet_normal.hpp"

//Unary
#include "gradient.hpp"
#include "divergence.hpp"

// Binary
#include "addition.hpp"
#include "inner_product.hpp"
#include "outer_product.hpp"
#include "colon_product.hpp"

namespace cfd
{

namespace forms
{

BilinearFormIntegral::cell_integral_tag dx;
BilinearFormIntegral::exterior_facet_integral_tag ds;
LinearForm::facet_normal_tag n;

forms::BilinearFormIntegral operator*(const forms::BilinearForm form, const forms::BilinearFormIntegral::cell_integral_tag region)
{
  return forms::BilinearFormIntegral(form, region);
}

forms::BilinearFormIntegral operator*(const forms::BilinearForm form, const forms::BilinearFormIntegral::exterior_facet_integral_tag region)
{
  return forms::BilinearFormIntegral(form, region);
}

forms::BilinearForm B(const forms::LinearForm trial, const forms::LinearForm test)
{
  return forms::BilinearForm(trial, test);
}

forms::LinearForm grad(const forms::LinearForm f)
{
  return forms::LinearForm(new forms::Gradient(f.getField()));
}

forms::LinearForm div(const forms::LinearForm f)
{
  return forms::LinearForm(new forms::Divergence(f.getField()));
}

forms::LinearForm operator+(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::Addition(l.getField(), r.getField()));
}

forms::BilinearFormIntegralSum operator+(const forms::BilinearFormIntegralSum l, const forms::BilinearFormIntegralSum r)
{
  forms::BilinearFormIntegralSum result(l);
  result.append(r);
  return result;
}

forms::LinearForm inner(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::InnerProduct(l.getField(), r.getField()));
}

forms::LinearForm outer(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::OuterProduct(l.getField(), r.getField()));
}

forms::LinearForm operator*(const forms::LinearForm l, const forms::LinearForm r)
{
  return outer(l, r);
}

forms::LinearForm colon(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::ColonProduct(l.getField(), r.getField()));
}

}

}

#endif
