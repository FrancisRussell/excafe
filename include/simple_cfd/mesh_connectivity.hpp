#ifndef SIMPLE_CFD_MESH_CONNECTIVITY_HPP
#define SIMPLE_CFD_MESH_CONNECTIVITY_HPP

#include<set>
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
  std::size_t numRelations(const std::size_t entity) const;
  std::size_t numRelations() const;
  void clear();
  std::vector<std::size_t> getIndices(const std::size_t entity) const;

  template<typename InputIterator>
  std::size_t addEntity(const InputIterator& indicesBegin, const InputIterator& indicesEnd)
  {
    const std::size_t index = offsets.size() - 1;
    addEntity(index, indicesBegin, indicesEnd);
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

    std::set<std::size_t> entities(indicesBegin, indicesEnd);
    entities.insert(indices.begin() + offsets[index], indices.begin() + offsets[index+1]); 
    const std::size_t numNewIndices = entities.size() + offsets[index] - offsets[index+1];

    // Now we make space for the new indices
    indices.insert(indices.begin() + offsets[index], numNewIndices, 0);
    std::copy(entities.begin(), entities.end(), indices.begin() + offsets[index]);

    // Now we have to increment all offsets after index by numNewIndices
    std::transform(offsets.begin() + index + 1, offsets.end(), offsets.begin() + index + 1, 
      boost::bind(std::plus<std::size_t>(), _1, numNewIndices));
  }

};

}

#endif
