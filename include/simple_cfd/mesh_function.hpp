#ifndef SIMPLE_CFD_MESH_FUNCTION_HPP
#define SIMPLE_CFD_MESH_FUNCTION_HPP

#include <cstddef>
#include <cassert>
#include <map>
#include <mesh_entity.hpp>

namespace cfd
{

template<typename T>
class MeshFunction
{
public:
  typedef T value_type;
  std::size_t dimension;
  value_type defaultValue;
  std::map<MeshEntity, value_type> relations;

  MeshFunction(const std::size_t _dimension, const value_type _defaultValue = value_type()) :
    dimension(_dimension), defaultValue(_defaultValue)
  {
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  void setValue(const MeshEntity& entity, const value_type value)
  {
    assert(entity.getDimension() == dimension);
    relations[entity] = value;
  }

  const value_type operator()(const MeshEntity& entity) const
  {
    assert(entity.getDimension() == dimension);
    const typename std::map<MeshEntity, value_type>::const_iterator entityIter(relations.find(entity));

    if (entityIter != relations.end())
    {
      return entityIter->second;
    }
    else
    {
      return defaultValue;
    }
  }
};

}

#endif
