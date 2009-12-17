#ifndef SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_EXPRESSION_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_ARRAY_EXPRESSION_HPP

#include <cstddef>
#include <set>
#include <boost/shared_ptr.hpp>

namespace cfd
{

namespace detail
{

class ArrayExpression
{
public:
  typedef boost::shared_ptr<ArrayExpression> expr_ptr;

  virtual std::size_t getTensorRank() const = 0;
  virtual std::size_t getTensorDimension() const = 0;
  virtual std::size_t numArrayIndices() const = 0;
  virtual std::size_t getArrayDimension(const std::size_t index) const = 0;
  virtual std::set<ArrayExpression*> getDependencies() const = 0;
  virtual ~ArrayExpression() {}
};

}

}
#endif
