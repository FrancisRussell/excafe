#ifndef SIMPLE_CFD_FORMS_TENSOR_LITERAL_HPP
#define SIMPLE_CFD_FORMS_TENSOR_LITERAL_HPP

#include "field.hpp"
#include "holders.hpp"
#include <simple_cfd/numeric/tensor.hpp>

namespace cfd
{

namespace forms
{

class TensorLiteral : public Field
{
private:
  TensorHolder literal;

public:
  template<std::size_t D>
  TensorLiteral(const Tensor<D>& tensor) : literal(tensor)
  {
  }

  std::size_t getRank() const
  {
    return literal.getRank();
  }

  std::size_t getDimension() const
  {
    return literal.getDimension();
  }
  
};

}

}

#endif
