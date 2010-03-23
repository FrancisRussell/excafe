#ifndef SIMPLE_CFD_CAPTURE_TENSOR_INDEX_VARIABLE_HPP
#define SIMPLE_CFD_CAPTURE_TENSOR_INDEX_VARIABLE_HPP

#include <cstddef>
#include <string>
#include <simple_cfd/exception.hpp>
#include "tensor_fwd.hpp"
#include "traits.hpp"

namespace cfd
{

namespace detail
{

template<typename T>
class TaggedIndexVariable
{
private:
  std::string name;
  std::size_t limit;

public:
  TaggedIndexVariable(const std::string& _name, const std::size_t _limit) :
    name(_name), limit(_limit)
  {
  }

  std::string getName() const
  {
    return name;
  }

  bool operator<(const TaggedIndexVariable& i) const
  {
    if (name != i.name)
    {
      return name < i.name;
    }
    else
    {
      if (limit != i.limit)
        CFD_EXCEPTION("Indices with identical names but different limits detected.");

      return false;
    }
  }

  bool operator==(const TaggedIndexVariable& i) const
  {
    if (name != i.name)
    {
      return false;
    }
    else
    {
      if (limit != i.limit)
        CFD_EXCEPTION("Indices with identical names but different limits detected.");

      return true;
    }
  }

  std::size_t getLimit() const
  {
    return limit;
  }
};

}

}
#endif
