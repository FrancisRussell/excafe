#include <set>
#include <algorithm>
#include <mesh_topology.hpp>
#include <general_cell.hpp>

namespace cfd
{

MeshTopology::MeshTopology(const GeneralCell& _cell) : cell(_cell), dimension(cell.getDimension()), 
  relations(numConnectivityRelations(dimension))
{
}

std::size_t MeshTopology::numConnectivityRelations(const std::size_t dimension)
{
  return (dimension+1)*(dimension+1);
}

void MeshTopology::setBaseConnectivity(const MeshConnectivity& c)
{
  MeshConnectivity* const baseConnectivity = getConnectivityObject(dimension, 0);
  assert(baseConnectivity->numEntities() == 0);
  *baseConnectivity = c;
}

std::size_t MeshTopology::numEntities(const std::size_t d)
{
  return getConnectivity(d, 0)->numEntities();
}

std::size_t MeshTopology::numRelations(const MeshEntity& entity, const std::size_t d)
{
  calculateConnectivity(entity.getDimension(), d);
  return getConnectivityObject(entity.getDimension(), d)->numRelations(entity.getIndex());
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
    performIntersection(d, dPrime, dPrimePrime);
  }
}

void MeshTopology::performTranspose(const std::size_t d, const std::size_t dPrime)
{
  MeshConnectivity* const newConnectivity = getConnectivityObject(d, dPrime);
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
  MeshConnectivity* const newConnectivity = getConnectivityObject(d, dPrime);
  assert(newConnectivity->numEntities() == 0);

  for(global_iterator dIter(global_begin(d)); dIter!=global_end(d); ++dIter)
  {
    const std::set<std::size_t> dVertexIndices(getIndices(*dIter, 0));

    for(local_iterator dPrimePrimeIter(local_begin(*dIter, dPrimePrime)); dPrimePrimeIter!=local_end(*dIter, dPrimePrime); ++dPrimePrimeIter) {
      for(local_iterator dPrimeIter(local_begin(*dPrimePrimeIter, dPrime)); dPrimeIter!=local_end(*dPrimePrimeIter, dPrime); ++dPrimeIter)
      {
        const std::set<std::size_t> dPrimeVertexIndices(getIndices(*dPrimeIter, 0));
        const bool dIncludesDPrime = std::includes(dVertexIndices.begin(), dVertexIndices.end(),
          dPrimeVertexIndices.begin(), dPrimeVertexIndices.end());

        if ((d == dPrime && dIter->getIndex() != dPrimeIter->getIndex()) || (d>dPrime && dIncludesDPrime)) 
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
  MeshConnectivity* const newConnectivity = getConnectivityObject(d, 0);
  assert(newConnectivity->numEntities() == 0);

  std::size_t k = 0;

  for(global_iterator cellIter(global_begin(dimension)); cellIter!=global_end(dimension); ++cellIter)
  {
    const std::set< std::set<std::size_t> > vi(cell.getIncidentVertices(*this, *cellIter, 0));
    for(local_iterator iCellIter(local_begin(*cellIter, dimension)); iCellIter!=local_end(*cellIter, dimension); ++iCellIter) 
    {
      if (iCellIter->getIndex() < cellIter->getIndex())
      {
        const std::set< std::set<std::size_t> > vj(cell.getIncidentVertices(*this, *iCellIter, 0));
        for(std::set< std::set<std::size_t> >::const_iterator viIter(vi.begin()); viIter!=vi.end(); ++viIter)
        {
          if (vj.find(*viIter) == vj.end())
          {
            newConnectivity->addEntity(k, viIter->begin(), viIter->end());
            ++k;
          }
        }
      }
    }
  }
}

std::size_t MeshTopology::getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const
{
  assert(0 <= d  && d <= dimension);
  assert(0 <= dPrime && dPrime <= dimension);
  return d * dimension + dPrime;
}

MeshConnectivity* MeshTopology::getConnectivityObject(const std::size_t d, const std::size_t dPrime)
{
  const std::size_t connectivityIndex = getConnectivityIndex(d, dPrime);
  return &relations[connectivityIndex];
}

MeshConnectivity* MeshTopology::getConnectivity(const std::size_t d, const std::size_t dPrime)
{
  calculateConnectivity(d, dPrime);
  return getConnectivityObject(d, dPrime);
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

MeshTopology::global_iterator MeshTopology::global_begin(const std::size_t d)
{
  return MeshEntityIteratorGlobal(this, d, 0);
}

MeshTopology::global_iterator MeshTopology::global_end(const std::size_t d)
{
  return MeshEntityIteratorGlobal(this, d, getConnectivity(d, 0)->numEntities());
}

MeshTopology::local_iterator MeshTopology::local_begin(const MeshEntity& entity, const std::size_t d)
{
  return MeshEntityIteratorLocal(this, entity, d, 0);
}

MeshTopology::local_iterator MeshTopology::local_end(const MeshEntity& entity, const std::size_t d)
{
  return MeshEntityIteratorLocal(this, entity, d, 
    getConnectivity(entity.getDimension(), d)->numRelations(entity.getIndex()));
}

}
