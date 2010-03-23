#ifndef SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_FWD_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_TENSOR_FWD_HPP

namespace cfd
{

namespace detail
{

template<typename T> class IndexExpression;
template<typename T> class TaggedIndexVariable;
template<typename T> class Index;

class ArraySize;
class TensorSize;
class TensorArrayTablePolynomial;
class IndexGenerator;
class TensorPlaceholder;
class TensorArrayPlaceholder;

// Tag classes
struct array_tag {};
struct tensor_tag {};
struct row_major_tag {};

// Typedefs for easier naming
typedef IndexExpression<array_tag>      ArrayIndexExpression;
typedef IndexExpression<tensor_tag>     TensorIndexExpression;
typedef TaggedIndexVariable<array_tag>  ArrayIndexVariable;
typedef TaggedIndexVariable<tensor_tag> TensorIndexVariable;
typedef Index<array_tag>                ArrayIndex;
typedef Index<tensor_tag>               TensorIndex;

}

}

#endif
