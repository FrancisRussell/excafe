#ifndef SIMPLE_CFD_FORMS_FORMS_HPP
#define SIMPLE_CFD_FORMS_FORMS_HPP

#include <cstddef>
#include <cassert>
#include <boost/shared_ptr.hpp>

#include "field.hpp"
#include "linear_form.hpp"

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

forms::LinearForm inner(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::InnerProduct(l.getField(), r.getField()));
}

forms::LinearForm outer(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::OuterProduct(l.getField(), r.getField()));
}

forms::LinearForm colon(const forms::LinearForm l, const forms::LinearForm r)
{
  return forms::LinearForm(new forms::ColonProduct(l.getField(), r.getField()));
}


}

#endif
