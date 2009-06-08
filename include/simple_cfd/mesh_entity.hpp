#ifndef SIMPLE_CFD_MESH_ENTITY_HPP
#define SIMPLE_CFD_MESH_ENTITY_HPP

namespace cfd
{

class MeshEntity
{
private:
  std::size_t dimension;
  std::size_t index;

public:
  MeshEntity(const std::size_t _dimension, const std::size_t _index) : dimension(_dimension), index(_index)
  {
  }

  std::size_t getDimension() const
  {
    return dimension;
  }

  std::size_t getIndex() const
  {
    return index;
  }

  MeshEntity& operator=(const MeshEntity& e)
  {
    dimension = e.dimension;
    index = e.index;
    return *this;
  }

  bool operator==(const MeshEntity& e) const
  {
    return dimension == e.dimension && index == e.index;
  }

  bool operator<(const MeshEntity& e) const
  {
    if (dimension < e.dimension) return true;
    if (dimension == e.dimension && index < e.index) return true;
    return false;
  }
};

}

#endif
