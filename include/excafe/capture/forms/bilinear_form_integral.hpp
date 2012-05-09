#ifndef EXCAFE_FORMS_BILINEAR_FORM_INTEGRAL_HPP
#define EXCAFE_FORMS_BILINEAR_FORM_INTEGRAL_HPP

#include <boost/variant.hpp>
#include <boost/static_assert.hpp>
#include "bilinear_form.hpp"

namespace excafe
{

namespace forms
{

class BilinearFormIntegral
{
public:
  class cell_integral_tag {};
  class exterior_facet_integral_tag {};
  typedef boost::variant<cell_integral_tag, exterior_facet_integral_tag> region_t;

private:
  BilinearForm form;
  region_t region;

public:
  template<typename integral_type>
  BilinearFormIntegral(const BilinearForm& _form, const integral_type& _region) : 
    form(_form), region(_region)
  {
  }

  BilinearForm getForm() const
  {
    return form;
  }

  region_t getRegion() const
  {
    return region;
  }
};

}

}

#endif
