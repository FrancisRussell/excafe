#ifndef SIMPLE_CFD_FORMS_BILINEAR_FORM_INTEGRAL_HPP
#define SIMPLE_CFD_FORMS_BILINEAR_FORM_INTEGRAL_HPP

#include "bilinear_form.hpp"

namespace cfd
{

namespace forms
{

class BilinearFormIntegral
{
public:
  enum Region
  {
    CELL,
    EXTERIOR_FACET
  };

private:
  BilinearForm form;
  Region region;

public:
  BilinearFormIntegral(const BilinearForm& _form, const Region _region) :
    form(_form), region(_region)
  {
    assert(region == CELL || region == EXTERIOR_FACET);
  }

  BilinearForm getForm() const
  {
    return form;
  }

  bool isCellIntegral() const
  {
    return region == CELL;
  }

  bool isExteriorFacetIntegral() const
  {
    return region == EXTERIOR_FACET;
  }
};

}

}

#endif
