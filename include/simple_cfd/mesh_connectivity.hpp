#ifndef SIMPLE_CFD_MESH_CONNECTIVITY_HPP
#define SIMPLE_CFD_MESH_CONNECTIVITY_HPP

#include<vector>
#include<cstddef>

namespace cfd
{

class MeshConnectivity
{
private:
  std::vector<std::size_t> indices;
  std::vector<std::size_t> offsets;

public:
  MeshConnectivity();
  std::size_t numEntities() const;
  std::size_t addEntity(const std::vector<std::size_t>& _indices);
  std::size_t entitySize(const std::size_t entity) const;
  void populateWithIndices(std::vector<std::size_t>& indices, const std::size_t entity) const;

  template<typename InputIterator>
  std::size_t addEntity(const InputIterator& indicesBegin, const InputIterator& indicesEnd)
  {
    const std::size_t index = offsets.size() - 1;
    indices.insert(indices.end(), indicesBegin, indicesEnd);
    offsets.push_back(indices.size());
    return index;
  }
};

}

#endif
