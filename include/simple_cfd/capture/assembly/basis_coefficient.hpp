#ifndef SIMPLE_CFD_CAPTURE_ASSEMBLY_BASIS_COEFFICIENT_HPP
#define SIMPLE_CFD_CAPTURE_ASSEMBLY_BASIS_COEFFICIENT_HPP

#include <cstddef>
#include <utility>
#include <simple_cfd/capture/fields/field.hpp>
#include <simple_cfd/capture/fields/discrete_field_expr.hpp>
#include "scalar_placeholder_operators.hpp"

namespace cfd
{

namespace detail
{

class BasisCoefficient : public ScalarPlaceholderOperators<BasisCoefficient>
{
private:
  DiscreteFieldExpr::expr_ptr fieldExpr;
  std::size_t index;

public:
  BasisCoefficient(const Field& field, const std::size_t _index) : 
    fieldExpr(field.getExpr()), index(_index)
  {
  }

  bool operator==(const BasisCoefficient& c) const
  {
    return fieldExpr == c.fieldExpr && index == c.index;
  }

  bool operator<(const BasisCoefficient& c) const
  {
    return std::make_pair(fieldExpr, index) <
      std::make_pair(c.fieldExpr, c.index);
  }

  DiscreteFieldExpr::expr_ptr getField() const
  {
    return fieldExpr;
  }

  std::size_t getIndex() const
  {
    return index;
  }
};

}

}

#endif
