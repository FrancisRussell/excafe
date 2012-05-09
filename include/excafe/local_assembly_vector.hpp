#ifndef EXCAFE_LOCAL_ASSEMBLY_VECTOR_HPP
#define EXCAFE_LOCAL_ASSEMBLY_VECTOR_HPP

#include <cstddef>
#include <set>
#include <boost/foreach.hpp>

namespace excafe
{

namespace detail
{

template<std::size_t D, typename T>
class LocalAssemblyVector
{
public:
  static const std::size_t dimension = D;
  typedef T value_type;
  typedef const FiniteElement<dimension> element_t;

private:
  std::set<const element_t*> elements;
  std::vector<value_type> values;

public:
  LocalAssemblyVector(const std::set<const element_t*>& _elements) : elements(_elements),
    values(getSize())
  {
  }

  std::size_t getSize() const
  {
    std::size_t size = 0;
    BOOST_FOREACH(const element_t* element, elements)
    {
      size += element->spaceDimension();
    }

    return size;
  }

  element_t& operator[](const std::size_t index)
  {
    assert(index < values.size());
    return values[index];
  }
};

}

}

#endif
