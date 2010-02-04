#ifndef SIMPLE_CFD_CAPTURE_ARRAY_FREE_TENSOR_ARRAY_HPP
#define SIMPLE_CFD_CAPTURE_ARRAY_FREE_TENSOR_ARRAY_HPP

#include <cstddef>

namespace cfd
{

namespace detail
{

class FreeTensorArray
{
private:
  std::size_t id;

public:
  FreeTensorArray(const std::size_t _id) : id(_id)
  {
  }

  std::size_t getID() const
  {
    return id;
  }

  bool operator==(const FreeTensorArray& a) const
  {
    return id == a.id;
  }

  bool operator<(const FreeTensorArray& a) const
  {
    return id < a.id;
  }
};

}

}

#endif
