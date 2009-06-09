#include <set>
#include <algorithm>
#include <mesh_topology.hpp>

namespace cfd
{

MeshTopology::MeshTopology(const std::size_t _dimension) : dimension(_dimension), 
  relations(numRelations(dimension))
{
}

std::size_t MeshTopology::numRelations(const std::size_t dimension)
{
  return (dimension+1)*(dimension+1);
}

void MeshTopology::calculateConnectivity(const std::size_t d, const std::size_t dPrime)
{
  const std::size_t connectivityIndex = getConnectivityIndex(d, dPrime);
  const MeshConnectivity& connectivity = relations[connectivityIndex];

  // Return if connectivity has already been calculated
  if (connectivity.numEntities() != 0)
  {
    return;
  }

  if (d < dPrime) // We want to calculate from the transpose
  {
    calculateConnectivity(dPrime, d);
    performTranspose(d, dPrime);
  }
  else if (dPrime == 0 && dimension > d && d > 0)
  {
    performBuild(d);
  }
  else
  {
    // We need to avoid the case d=d'=d''. We choose d''=0 and have d''=dimension
    // in the special case where d=d'=0.
    const std::size_t dPrimePrime = (d == 0 && dPrime == 0) ? dimension : 0;
    calculateConnectivity(d, dPrimePrime);
    calculateConnectivity(dPrimePrime, dPrime);
    performIntersection(d, dPrime, dPrimePrime);
  }
}

void MeshTopology::performTranspose(const std::size_t d, const std::size_t dPrime)
{
  MeshConnectivity* const newConnectivity = getConnectivity(d, dPrime);
  assert(newConnectivity->numEntities() == 0);

  for(global_iterator dPrimeIter(global_begin(dPrime)); dPrimeIter!=global_end(dPrime); ++dPrimeIter)
  {
    for(local_iterator dIter(local_begin(*dPrimeIter, d)); dIter!=local_end(*dPrimeIter, d); ++dIter)
    {
      const std::size_t j = dPrimeIter->getIndex();
      newConnectivity->addEntity(dIter->getIndex(), &j, &j+1);
    }
  }
}

void MeshTopology::performIntersection(const std::size_t d, const std::size_t dPrime, const std::size_t dPrimePrime)
{
  calculateConnectivity(dPrime, 0);
  calculateConnectivity(d, 0);

  MeshConnectivity* const newConnectivity = getConnectivity(d, dPrime);
  assert(newConnectivity->numEntities() == 0);

  for(global_iterator dIter(global_begin(d)); dIter!=global_end(d); ++dIter)
  {
    const std::set<std::size_t> dVertexIndices(getIndices(*dIter, 0));

    for(local_iterator dPrimePrimeIter(local_begin(*dIter, dPrimePrime)); dPrimePrimeIter!=local_end(*dIter, dPrimePrime); ++dPrimePrimeIter)
    {
      for(local_iterator dPrimeIter(local_begin(*dPrimePrimeIter, dPrime)); dPrimeIter!=local_end(*dPrimePrimeIter, dPrime); ++dPrimeIter)
      {
        const std::set<std::size_t> dPrimeVertexIndices(getIndices(*dPrimeIter, 0));
        const bool dPrimeContainsD = std::includes(dPrimeVertexIndices.begin(), dPrimeVertexIndices.end(), dVertexIndices.begin(), dVertexIndices.end());

        if ((d == dPrime && dIter->getIndex() != dPrimeIter->getIndex()) || (d>dPrime && dPrimeContainsD)) 
        {
           const std::size_t j = dPrimeIter->getIndex();
           newConnectivity->addEntity(dIter->getIndex(), &j, &j+1);
        }
      }
    }
  }
}

void MeshTopology::performBuild(const std::size_t d)
{
  calculateConnectivity(dimension, dimension);
  std::size_t k = 0;

  for(global_iterator dIter(global_begin(dimension)); dIter!=global_end(dimension); ++dIter)
  {
  }
}

std::size_t MeshTopology::getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const
{
  assert(0 <= d  && d <= dimension);
  assert(0 <= dPrime && dPrime <= dimension);
  return d * dimension + dPrime;
}

MeshConnectivity* MeshTopology::getConnectivity(const std::size_t d, const std::size_t dPrime)
{
  const std::size_t connectivityIndex = getConnectivityIndex(d, dPrime);
  return &relations[connectivityIndex];
}

std::set<std::size_t> MeshTopology::getIndices(const MeshEntity& entity, const std::size_t d)
{
  std::set<std::size_t> indices;

  for(local_iterator dIter(local_begin(entity, d)); dIter!=local_end(entity, d); ++dIter)
  {
    indices.insert(dIter->getIndex());
  }

  return indices;
}

}
