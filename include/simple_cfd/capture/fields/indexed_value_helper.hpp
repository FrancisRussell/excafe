#ifndef SIMPLE_CFD_CAPTURE_FIELDS_INDEXED_VALUE_HELPER_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_INDEXED_VALUE_HELPER_HPP

namespace cfd
{

namespace detail
{

template<typename discrete_object_type>
class IndexedValueHelper
{
private:
  typedef IndexableValue<discrete_object_type> parent_t
  parent_t::value_ptr parent;

public:
  IndexedValueHelper(const parent_t::value_ptr& _parent)
  {
  }

// Casting operator to actual discrete type
// assignment operator that refers to parent
};

}

};

#endif
