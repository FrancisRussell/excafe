#ifndef SIMPLE_CFD_MESH_FUNCTION_HPP
#define SIMPLE_CFD_MESH_FUNCTION_HPP

#include<cstddef>
#include<cassert>
#include<map>
#include<mesh_entity.hpp>

namespace cfd
{

template<typename T>
class MeshFunction
{
public:
  typedef T value_type;
  std::size_t dimension;
  std::map<MeshEntity, value_type> relations;

  MeshFunction(const std::size_t _dimension) : dimension(_dimension)
  {
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  value_type& operator()(const MeshEntity& entity)
  {
    assert(entity.getDimension() == dimension);
    return relations[entity];
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
      return value_type();
    }
  }
};

}

#endif
