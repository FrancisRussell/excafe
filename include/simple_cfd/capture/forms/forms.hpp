#ifndef SIMPLE_CFD_FORMS_FORMS_HPP
#define SIMPLE_CFD_FORMS_FORMS_HPP

#include <cstddef>
#include <cassert>
#include <boost/shared_ptr.hpp>

#include "field_expr.hpp"
#include "linear_form.hpp"
#include "bilinear_form.hpp"
#include "bilinear_form_integral.hpp"
#include "bilinear_form_integral_sum.hpp"
#include "facet_normal.hpp"

//Unary
#include "field_gradient.hpp"
#include "field_divergence.hpp"

// Binary
#include "field_addition.hpp"
#include "field_inner_product.hpp"
#include "field_outer_product.hpp"
#include "field_colon_product.hpp"

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
  return forms::LinearForm(new detail::FieldGradient(f.getField()));
}

forms::LinearForm div(const forms::LinearForm f)
{
  return forms::LinearForm(new detail::FieldDivergence(f.getField()));
}

forms::LinearForm operator+(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new detail::FieldAddition(l.getField(), r.getField()));
}

forms::BilinearFormIntegralSum operator+(const forms::BilinearFormIntegralSum l, const forms::BilinearFormIntegralSum r)
{
  forms::BilinearFormIntegralSum result(l);
  result.append(r);
  return result;
}

forms::LinearForm inner(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new detail::FieldInnerProduct(l.getField(), r.getField()));
}

forms::LinearForm outer(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new detail::FieldOuterProduct(l.getField(), r.getField()));
}

forms::LinearForm operator*(const forms::LinearForm l, const forms::LinearForm r)
{
  return outer(l, r);
}

forms::LinearForm colon(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new detail::FieldColonProduct(l.getField(), r.getField()));
}

}

}

#endif
