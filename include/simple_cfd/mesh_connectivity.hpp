#ifndef SIMPLE_CFD_MESH_CONNECTIVITY_HPP
#define SIMPLE_CFD_MESH_CONNECTIVITY_HPP

#include<vector>
#include<algorithm>
#include<cstddef>
#include<functional>
#include<iterator>
#include<boost/bind.hpp>

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
  void clear();
  void populateWithIndices(std::vector<std::size_t>& indices, const std::size_t entity) const;

  template<typename InputIterator>
  std::size_t addEntity(const InputIterator& indicesBegin, const InputIterator& indicesEnd)
  {
    const std::size_t index = offsets.size() - 1;
    indices.insert(indices.end(), indicesBegin, indicesEnd);
    offsets.push_back(indices.size());
    return index;
  }

  template<typename InputIterator>
  void addEntity(const std::size_t index, const InputIterator& indicesBegin, const InputIterator& indicesEnd)
  {
    // First we need to resize the array of offsets if index is larger than the current largest entity
    // Be careful with unsigned types here!
    if (numEntities() <= index)
    {
      const std::size_t newOffsetsCount = index - numEntities() + 1;
      std::fill_n(std::back_inserter(offsets), newOffsetsCount, indices.size());
    }

    // Now we insert the new indices
    // We use index+1 so we insert new indices after existing indices for that entity
    const std::size_t numNewIndices = std::distance(indicesBegin, indicesEnd);
    indices.insert(indices.begin() + offsets[index+1], indicesBegin, indicesEnd);

    // Now we have to increment all offsets after index by numNewIndices
    std::transform(offsets.begin() + index + 1, offsets.end(), offsets.begin() + index + 1, 
      boost::bind(std::plus<std::size_t>(), _1, numNewIndices));
  }

};

}

#endif
