#ifndef SIMPLE_CFD_CAPTURE_TENSOR_INDEX_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_INDEX_HPP

#include <cstddef>
#include <vector>
#include "tensor_fwd.hpp"
#include "traits.hpp"

namespace cfd
{

namespace detail
{

template<typename T>
class Index
{
private:
  typedef T object_tag;
  typedef typename IndexProperties<object_tag>::size_type size_type;
  typedef typename IndexProperties<object_tag>::index_expression_t index_expression_t; 

  size_type size;
  std::vector<index_expression_t> indices;

public:
  typedef typename index_expression_t::constant_t constant_t;
  typedef typename index_expression_t::index_variable_t index_variable_t;

  Index(const size_type& _size) : size(_size), indices(size.numIndices())
  {
  }

  index_expression_t& operator[](const std::size_t index)
  {
    return indices[index];
  }

  const index_expression_t& operator[](const std::size_t index) const
  {
    return indices[index];
  }
};

}

}

#endif
