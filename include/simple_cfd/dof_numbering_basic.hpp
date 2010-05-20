#ifndef SIMPLE_CFD_DOF_NUMBERING_BASIC_HPP
#define SIMPLE_CFD_DOF_NUMBERING_BASIC_HPP

#include <cstddef>
#include <boost/array.hpp>
#include <boost/scoped_ptr.hpp>
#include <numeric>
#include <utility>
#include <cassert>
#include <map>
#include "mesh.hpp"
#include "general_cell.hpp"
#include "mesh_entity.hpp"
#include "dof_association.hpp"

namespace cfd
{

template<std::size_t D>
class DofNumberingBasic
{
private:
  static const std::size_t dimension = D;

  boost::scoped_ptr< GeneralCell<dimension> > cell;
  const boost::array<std::size_t, dimension+1> dofsPerEntity;
  const boost::array<std::size_t, dimension+1> numCellEntities;
  const std::size_t tensorSize;

  static boost::array<std::size_t, dimension+1> buildNumCellEntities(const GeneralCell<dimension>& cell)
  {
    boost::array<std::size_t, dimension+1> numCellEntities;

    for(std::size_t i=0; i<=dimension; ++i)
      numCellEntities[i] = cell.numEntities(i);

    return numCellEntities;
  }

  std::size_t numDofsPerValue() const
  {
    return std::inner_product(dofsPerEntity.begin(), dofsPerEntity.end(), numCellEntities.begin(), 0);
  }

public:
  DofNumberingBasic(const GeneralCell<dimension>& _cell, const boost::array<std::size_t, dimension+1>& _dofsPerEntity,
    const std::size_t _tensorSize) : cell(_cell.cloneGeneralCell()), dofsPerEntity(_dofsPerEntity),
    numCellEntities(buildNumCellEntities(*cell)), tensorSize(_tensorSize)
  {
  }

  DofNumberingBasic(const DofNumberingBasic& d) : cell(d.cell->cloneGeneralCell()),
    dofsPerEntity(d.dofsPerEntity), numCellEntities(d.numCellEntities), tensorSize(d.tensorSize)
  {
  }

  std::size_t getTensorIndex(const std::size_t dof) const
  {
    assert(dof < numDofs());
    const std::size_t dofsPerValue = numDofsPerValue();
    return dof / dofsPerValue;
  }

  DofAssociation getLocalAssociation(const std::size_t dof) const
  {
    assert(dof < numDofs());
    const std::size_t dofsPerValue = numDofsPerValue();
    const std::size_t scalarDof = dof % dofsPerValue;

    std::size_t entityDof = scalarDof;
    for(std::size_t i=0; i<=dimension; ++i)
    {
      if (entityDof < dofsPerEntity[i] * numCellEntities[i])
      {
        const MeshEntity localEntity(i, entityDof / dofsPerEntity[i]);
        const std::size_t localIndex(entityDof % dofsPerEntity[i]);

        return DofAssociation(localEntity, localIndex);
      }
      else
      {
        entityDof -= dofsPerEntity[i] * numCellEntities[i];
      }
    }

    assert(false && "Index out of range");
    return DofAssociation(MeshEntity(0, 0), 0);
  }

  std::size_t numDofs() const
  {
    return numDofsPerValue() * tensorSize;
  }
};

}

#endif
