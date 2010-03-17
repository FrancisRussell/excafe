#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_FWD_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_FWD_HPP

namespace cfd
{

namespace detail
{

template<typename T> class IndexExpression;
template<typename T> class TaggedIndexVariable;

// Tag classes
struct array_tag {};
struct tensor_tag {};

// Typedefs for easier naming
typedef IndexExpression<array_tag>      ArrayIndexExpression;
typedef IndexExpression<tensor_tag>     TensorIndexExpression;
typedef TaggedIndexVariable<array_tag>  ArrayIndexVariable;
typedef TaggedIndexVariable<tensor_tag> TensorIndexVariable;

}

}

#endif
