#ifndef SIMPLE_CFD_CAPTURE_FIELDS_BOUNDARY_CONDITION_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_BOUNDARY_CONDITION_HPP

#include <cstddef>

namespace cfd
{

class BoundaryCondition
{
private:
  std::size_t index;

public:
  BoundaryCondition()
  {
  }

  BoundaryCondition(const std::size_t _index) : index(_index)
  {
  }

  std::size_t getIndex() const
  {
    return index;
  }

  BoundaryCondition& operator=(const BoundaryCondition& b)
  {
    index = b.index;
    return *this;
  }
};

}

#endif
