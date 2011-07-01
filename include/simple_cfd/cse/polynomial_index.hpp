#ifndef SIMPLE_CFD_CSE_POLYNOMIAL_INDEX_HPP
#define SIMPLE_CFD_CSE_POLYNOMIAL_INDEX_HPP

#include <ostream>
#include <limits>
#include <boost/operators.hpp>

namespace cfd
{

namespace cse
{

class PolynomialIndex : boost::totally_ordered<PolynomialIndex>
{
private:
  std::size_t index;

public:
  PolynomialIndex() : index(std::numeric_limits<std::size_t>::max())
  {
  }

  explicit PolynomialIndex(const std::size_t _index) : index(_index)
  {
  }

  std::size_t getIndex() const
  {
    return index;
  }

  bool operator==(const PolynomialIndex& i) const
  {
    return index == i.index;
  }

  bool operator<(const PolynomialIndex& i) const
  {
    return index < i.index;
  }

  void write(std::ostream& o) const
  {
    o << "var(" << index << ")";
  }
};

std::ostream& operator<<(std::ostream& o, const PolynomialIndex& i);

}

}

#endif
