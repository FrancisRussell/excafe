#ifndef EXCAFE_CAPTURE_ASSEMBLY_BASIS_COEFFICIENT_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_BASIS_COEFFICIENT_HPP

#include <cstddef>
#include <utility>
#include <ostream>
#include <excafe/capture/fields/field.hpp>
#include <excafe/capture/fields/discrete_field_expr.hpp>

namespace excafe
{

namespace detail
{

class BasisCoefficient
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

  void write(std::ostream& o) const
  {
    o << "field[" << fieldExpr << "][" << index << "]";
  }
};

std::ostream& operator<<(std::ostream& o, const BasisCoefficient& c);

}

}

#endif
