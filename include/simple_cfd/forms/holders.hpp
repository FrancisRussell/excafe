#ifndef SIMPLE_CFD_FORMS_HOLDERS_HPP
#define SIMPLE_CFD_FORMS_HOLDERS_HPP

#include <cstddef>
#include <cassert>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <simple_cfd/mesh.hpp>
#include <simple_cfd/finite_element.hpp>
#include <simple_cfd/dof_map.hpp>
#include <simple_cfd/numeric/tensor.hpp>

namespace cfd
{

namespace forms
{

class FiniteElementHolder
{
private:
  boost::any elementPtr;
  std::size_t dimension;
  std::size_t rank;

public:
  template<std::size_t D>
  FiniteElementHolder(const FiniteElement<D>& _element) : elementPtr(&_element), 
    dimension(_element.getDimension()), rank(_element.getRank())
  {
  }

  boost::any getElementPtr() const
  {
    return elementPtr;
  }

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }
};

class FEVectorHolder
{
private:
  boost::any vectorPtr;
  std::size_t rank;
  std::size_t dimension;

public:
  template<std::size_t D>
  FEVectorHolder(const FEVector<D>& vector) : vectorPtr(&vector), rank(vector.getRank()),
    dimension(vector.getDimension())
  {
  }

  boost::any getVectorPtr() const
  {
    return vectorPtr;
  }

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }
};

class TensorHolder
{
private:
  boost::any tensor;
  std::size_t dimension;
  std::size_t rank;

public:
  template<std::size_t D>
  TensorHolder(const Tensor<D>& t) : tensor(t), 
  dimension(t.getDimension()), rank(t.getRank())
  {
  }

  std::size_t getRank() const
  {
    return rank;
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  boost::any getTensor() const
  {
    return tensor;
  }
};

}

}

#endif
