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

  calculateConnectivity(dimension, dimension);
}

void MeshTopology::performIntersection(const std::size_t d, const std::size_t dPrime, const std::size_t dPrimePrime)
{
}

void MeshTopology::performBuild(const std::size_t d)
{
}

std::size_t MeshTopology::getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const
{
  assert(0 <= d  && d <= dimension);
  assert(0 <= dPrime && dPrime <= dimension);
  return d * dimension + dPrime;
}

}
