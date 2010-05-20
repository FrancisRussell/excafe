#ifndef SIMPLE_CFD_DOF_ASSOCIATION_HPP
#define SIMPLE_CFD_DOF_ASSOCIATION_HPP

namespace cfd
{

class DofAssociation
{
private:
  MeshEntity entity;
  std::size_t index;

public:
  DofAssociation(const MeshEntity& _entity, const std::size_t _index) : entity(_entity), index(_index)
  {
  }

  MeshEntity getEntity() const
  {
    return entity;
  }

  std::size_t getEntityDimension() const
  {
    return entity.getDimension();
  }

  std::size_t getEntityIndex() const
  {
    return entity.getIndex();
  }

  std::size_t getIndex() const
  {
    return index;
  }

  bool operator==(const DofAssociation& a)
  {
    return entity == a.entity && index == a.index;
  }

  bool operator<(const DofAssociation& a)
  {
    if (entity < a.entity) return true;
    if (entity == a.entity && index < a.index) return true;
    return false;
  }

  DofAssociation& operator=(const DofAssociation& a)
  {
    entity = a.entity;
    index = a.index;
    return *this;
  }
};

}

#endif
