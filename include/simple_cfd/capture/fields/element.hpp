#ifndef SIMPLE_CFD_CAPTURE_FIELDS_ELEMENT_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_ELEMENT_HPP

#include <cstddef>

namespace cfd
{

class Element
{
private:
  std::size_t index;

public:
  Element()
  {
  }

  Element(const std::size_t index)
  {
  }

  std::size_t getIndex() const
  {
    return index;
  }
};

}

#endif
