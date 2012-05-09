#ifndef EXCAFE_BOUNDARY_CONDITION_TRIVIAL_HPP
#define EXCAFE_BOUNDARY_CONDITION_TRIVIAL_HPP

#include "boundary_condition3.hpp"

namespace excafe
{

template<std::size_t D>
class BoundaryConditionTrivial : public BoundaryCondition3<D>
{
private:
  static const std::size_t dimension = D;
  const int label;
  const Tensor<dimension> value;

public:
  BoundaryConditionTrivial(const int _label, const Tensor<dimension>& _value) : label(_label), value(_value)
  {
  }

  virtual std::size_t getRank() const
  {
    return value.getRank();
  }

  virtual bool applies(const int l) const
  {
    return label == l;
  }

  virtual Tensor<dimension> getValue(const vertex<dimension>& location, const int label) const
  {
    return value;
  }
};

}

#endif
