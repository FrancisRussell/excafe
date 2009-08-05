#ifndef SIMPLE_CFD_FORM_BASIS_FINDER_HPP
#define SIMPLE_CFD_FORM_BASIS_FUNDER_HPP

#include <cstddef>
#include <set>
#include <boost/any.hpp>
#include "field_visitor.hpp"
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/discrete_field.hpp>

namespace cfd
{

namespace forms
{

template<std::size_t D> 
class BasisFinder : public FieldVisitor
{
private:
  static const std::size_t dimension = D;
  typedef FiniteElement<dimension> finite_element_t;
  const finite_element_t* basis;

  void handle(const FiniteElementHolder& holder)
  {
    const finite_element_t* const elementPtr = boost::any_cast<const finite_element_t*>(holder.getElementPtr());
    assert(elementPtr != NULL);
    handle(elementPtr);
  }

  void handle(const DiscreteFieldHolder& holder)
  {
    const DiscreteField<dimension>* const vectorPtr = boost::any_cast<const DiscreteField<dimension>*>(holder.getVectorPtr());
    assert(vectorPtr != NULL);

    const std::set<const finite_element_t*> elements(vectorPtr->getRowMappings().getFiniteElements());
    assert(elements.size() == 1);

    handle(*elements.begin());
  }

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
  BasisFinder() : basis(NULL)
  {
  }

  const finite_element_t* getBasis() const
  {
    return basis;
  }

  virtual void enter(Addition& addition)
  {
  }

  virtual void exit(Addition& addition)
  {
  }

  virtual void enter(InnerProduct& addition)
  {
  }

  virtual void exit(InnerProduct& addition)
  {
  }

  virtual void enter(OuterProduct& addition)
  {
  }

  virtual void exit(OuterProduct& addition)
  {
  }

  virtual void enter(ColonProduct& addition)
  {
  }

  virtual void exit(ColonProduct& addition)
  {
  }

  virtual void enter(Gradient& gradient)
  {
  }

  virtual void exit(Gradient& gradient)
  {
  }

  virtual void enter(Divergence& gradient)
  {
  }

  virtual void exit(Divergence& gradient)
  {
  }

  // Terminals
  virtual void visit(FacetNormal& normal)
  {
  }

  virtual void visit(BasisField& basis)
  {
    handle(basis.getElement());
  }

  virtual void visit(DiscreteFieldReference& d)
  {
    handle(d.getVector());
  }

  virtual void visit(TensorLiteral& tensor)
  {
    // We don't reference any fields
  }
};

}

}

#endif
