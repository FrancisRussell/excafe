#include <mesh_topology.hpp>
#include <general_cell.hpp>
#include <iterator>
#include <set>
#include <vector>
#include <algorithm>

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
  if (d>0)
  {
    // This will use build (except when d==dimension)
    return getConnectivity(d, 0)->numEntities();
  }
  else
  {
    // This will use transpose
    return getConnectivity(0, dimension)->numEntities();
  }
}

std::size_t MeshTopology::numRelations(const MeshEntity& entity, const std::size_t d)
{
  return getConnectivity(entity.getDimension(), d)->numRelations(entity.getIndex());
}

std::size_t MeshTopology::numRelations(const std::size_t d, const std::size_t dPrime)
{
  return getConnectivity(d, dPrime)->numRelations();
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
  else if (dPrime == 0 && d == dimension)
  {
    // We need D->0 to already exist!
    assert(false);
  }
  else if (dPrime == 0 && d > 0)
  {
    performBuild(d);
  }
  else if (dPrime == 0 && d == 0)
  {
    // We need this because 0->0 intersection via D is invalid
    performBuildZeroToZero();
  }
  else
  {
    performIntersection(d, dPrime, 0);
  }
}

void MeshTopology::performTranspose(const std::size_t d, const std::size_t dPrime)
{
  assert(d < dPrime);

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

bool MeshTopology::contains(const MeshEntity& m1, const MeshEntity& m2)
{
  if (m2.getDimension() > m1.getDimension())
  {
    return false;
  }
  else if (m1.getDimension() == m2.getDimension())
  {
    return m1.getIndex() == m2.getIndex();
  }
  else
  {
    std::set<std::size_t> m1VertexIndices;
    outputIndices(m1, 0, std::inserter(m1VertexIndices, m1VertexIndices.begin()));

    std::set<std::size_t> m2VertexIndices;
    outputIndices(m2, 0, std::inserter(m2VertexIndices, m2VertexIndices.begin()));

    const bool m1_includes_m2 = std::includes(m1VertexIndices.begin(), m1VertexIndices.end(),
                                m2VertexIndices.begin(), m2VertexIndices.end());
    return m1_includes_m2;
  }
}

void MeshTopology::performIntersection(const std::size_t d, const std::size_t dPrime, const std::size_t dPrimePrime)
{
  assert(d >= dPrime);
  MeshConnectivity* const newConnectivity = getConnectivityObject(d, dPrime);
  assert(newConnectivity->numEntities() == 0);

  for(global_iterator dIter(global_begin(d)); dIter!=global_end(d); ++dIter)
  {
    for(local_iterator dPrimePrimeIter(local_begin(*dIter, dPrimePrime)); dPrimePrimeIter!=local_end(*dIter, dPrimePrime); ++dPrimePrimeIter)
    {
      for(local_iterator dPrimeIter(local_begin(*dPrimePrimeIter, dPrime)); dPrimeIter!=local_end(*dPrimePrimeIter, dPrime); ++dPrimeIter)
      {
        if ((d == dPrime && dIter->getIndex() != dPrimeIter->getIndex()) || (d>dPrime && contains(*dIter, *dPrimeIter))) 
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
    const std::vector< std::set<std::size_t> > vi(cell.getIncidentVertices(*this, *cellIter, d));
    std::set< std::set<std::size_t> > seen;

    for(local_iterator iCellIter(local_begin(*cellIter, dimension)); iCellIter!=local_end(*cellIter, dimension); ++iCellIter)
    {
      if (iCellIter->getIndex() < cellIter->getIndex())
      {
        const std::vector< std::set<std::size_t> > vj(cell.getIncidentVertices(*this, *iCellIter, d));
        seen.insert(vj.begin(), vj.end());
      }
    }

    for(std::vector< std::set<std::size_t> >::const_iterator viIter(vi.begin()); viIter!=vi.end(); ++viIter)
    {
      if (seen.find(*viIter) == seen.end())
      {
        newConnectivity->addEntity(k, viIter->begin(), viIter->end());
        ++k;
      }
    }
  }
}

void MeshTopology::performBuildZeroToZero()
{
  MeshConnectivity* const newConnectivity = getConnectivityObject(0, 0);
  assert(newConnectivity->numEntities() == 0);

  for(global_iterator cellIter(global_begin(dimension)); cellIter!=global_end(dimension); ++cellIter)
  {
    for(local_iterator vIter(local_begin(*cellIter, 0)); vIter!=local_end(*cellIter, 0); ++vIter)
    {
      const std::size_t vid = vIter->getIndex();
      newConnectivity->addEntity(vid, &vid, &vid+1);
    }
  }
}


std::size_t MeshTopology::getConnectivityIndex(const std::size_t d, const std::size_t dPrime) const
{
  assert(0 <= d  && d <= dimension);
  assert(0 <= dPrime && dPrime <= dimension);
  return d * (dimension+1) + dPrime;
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

std::vector<std::size_t> MeshTopology::getIndices(const MeshEntity& entity, const std::size_t d)
{
  return getConnectivity(entity.getDimension(), d)->getIndices(entity.getIndex());
}

MeshTopology::global_iterator MeshTopology::global_begin(const std::size_t d)
{
  return MeshEntityIteratorGlobal(this, d, 0);
}

MeshTopology::global_iterator MeshTopology::global_end(const std::size_t d)
{
  return MeshEntityIteratorGlobal(this, d, numEntities(d));
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
