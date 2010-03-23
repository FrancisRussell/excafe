#ifndef SIMPLE_CFD_CAPTURE_TENSOR_INDEX_EXPRESSION_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_INDEX_EXPRESSION_HPP

#include <cstddef>
#include <boost/variant.hpp>
#include <simple_cfd/exception.hpp>
#include "tensor_fwd.hpp"

namespace cfd
{

namespace detail
{

template<typename T>
class IndexExpression
{
private:
  typedef T object_tag;

public:
  typedef std::size_t constant_t;
  typedef typename IndexProperties<object_tag>::index_variable_t index_variable_t;

private:
  typedef boost::variant<constant_t, index_variable_t> variant_t;
  variant_t value;

public:
  IndexExpression() : value(0)
  {
  }

  IndexExpression(const constant_t c) : value(c)
  {
  }

  IndexExpression(const index_variable_t i) : value(i)
  {
  }

  bool isConstant() const
  {
    return boost::get<constant_t>(&value) != NULL;
  }

  bool isVariable() const
  {
    return !isConstant();
  }

  constant_t toConstant() const
  {
    if (!isConstant())
      CFD_EXCEPTION("Attemped to convert an index variable valued index expression to a constant.");

    return boost::get<constant_t>(value);
  }

  index_variable_t toIndexVariable() const
  {
    if (!isVariable())
      CFD_EXCEPTION("Attemped to convert a constant valued index expression to an index variable.");

    return boost::get<index_variable_t>(value);
  }

  bool operator==(const IndexExpression& e) const
  {
    return value == e.value;
  }

  bool operator<(const IndexExpression& e) const
  {
    return value < e.value;
  }
};

}

}

#endif
