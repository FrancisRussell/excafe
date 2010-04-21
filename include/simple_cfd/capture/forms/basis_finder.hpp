#ifndef SIMPLE_CFD_FORM_BASIS_FINDER_HPP
#define SIMPLE_CFD_FORM_BASIS_FINDER_HPP

#include <cstddef>
#include <cassert>
#include <set>
#include "field_visitor.hpp"
#include "field_addition.hpp"
#include "field_inner_product.hpp"
#include "field_outer_product.hpp"
#include "field_colon_product.hpp"
#include "field_gradient.hpp"
#include "field_divergence.hpp"
#include "facet_normal.hpp"
#include "field_basis.hpp"
#include "field_discrete_reference.hpp"
#include "field_scalar.hpp"
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/discrete_field.hpp>
#include <simple_cfd/capture/evaluation/evaluation_fwd.hpp>
#include <simple_cfd/capture/capture_fwd.hpp>

namespace cfd
{

namespace forms
{

template<std::size_t D> 
class BasisFinder : public detail::FieldVisitor
{
private:
  static const std::size_t dimension = D;
  typedef FiniteElement<dimension> finite_element_t;

  Scenario<dimension>& scenario;
  const finite_element_t* basis;

  void handle(const finite_element_t* const element)
  {
    assert(element != NULL);

    if (basis == NULL || basis == element)
    {
      basis = element;
    }
    else
    {
      assert(false && "Field appears to have more than one basis!");
    }
  }

public:
  BasisFinder(Scenario<dimension>& _scenario) : scenario(_scenario),
   basis(NULL)
  {
  }

  const finite_element_t* getBasis() const
  {
    assert(basis != NULL);
    return basis;
  }

  virtual void enter(detail::FieldAddition& addition)
  {
  }

  virtual void exit(detail::FieldAddition& addition)
  {
  }

  virtual void enter(detail::FieldInnerProduct& addition)
  {
  }

  virtual void exit(detail::FieldInnerProduct& addition)
  {
  }

  virtual void enter(detail::FieldOuterProduct& addition)
  {
  }

  virtual void exit(detail::FieldOuterProduct& addition)
  {
  }

  virtual void enter(detail::FieldColonProduct& addition)
  {
  }

  virtual void exit(detail::FieldColonProduct& addition)
  {
  }

  virtual void enter(detail::FieldGradient& gradient)
  {
  }

  virtual void exit(detail::FieldGradient& gradient)
  {
  }

  virtual void enter(detail::FieldDivergence& gradient)
  {
  }

  virtual void exit(detail::FieldDivergence& gradient)
  {
  }

  // Terminals
  virtual void visit(detail::FacetNormal& normal)
  {
  }

  virtual void visit(detail::FieldBasis& basis)
  {
    handle(&scenario.getElement(basis.getElement()));
  }

  virtual void visit(detail::FieldDiscreteReference& d)
  {
  }

  virtual void visit(detail::FieldScalar& s)
  {
  }
};

}

}

#endif
