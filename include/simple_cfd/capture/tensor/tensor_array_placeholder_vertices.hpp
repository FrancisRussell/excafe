#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_VERTICES_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_ARRAY_PLACEHOLDER_VERTICES_HPP

namespace cfd
{

namespace detail
{

class TensorArrayPlaceholderVertices : public TensorArrayPlaceholder
{
public:
  ArrayIndexVariable vertexIndex;
public:
  TensorArrayPlaceholderVertices(const ArrayIndexVariable& _vertexIndex, const std::size_t dimension) :
    TensorArrayPlaceholder(TensorSize(1, dimension)), vertexIndex(_vertexIndex)
  {
  }

};

}

}

#endif
