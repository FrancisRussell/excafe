#ifndef SIMPLE_CFD_NUMERIC_FWD_HPP
#define SIMPLE_CFD_NUMERIC_FWD_HPP

namespace cfd
{

namespace detail
{

// Tags for object types
struct tensor_tag {};

// General index type
template<typename T> class Index;

}

// Polynomial types
template<typename V, typename T> class Monomial;
template<typename V> class Polynomial;
template<typename V> class PolynomialFraction;
template<typename V> class OptimisedPolynomial;
template<typename V> class OptimisedPolynomialFraction;
template<typename V> class GinacExpression;
template<typename V> class NumericExpression;
template<typename V> class NumericExpressionVisitor;

template<std::size_t D, typename T> class Tensor;
class TensorSize;
typedef detail::Index<detail::tensor_tag> TensorIndex;

}

#endif
