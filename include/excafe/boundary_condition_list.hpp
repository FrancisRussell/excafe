#ifndef EXCAFE_BOUNDARY_CONDITION_LIST_HPP
#define EXCAFE_BOUNDARY_CONDITION_LIST_HPP

#include <cstddef>
#include <cassert>
#include <vector>
#include <boost/shared_ptr.hpp>
#include "boundary_condition3.hpp"
#include "vertex.hpp"

namespace excafe
{

template<std::size_t D>
class BoundaryConditionList : public BoundaryCondition3<D>
{
private:
  static const std::size_t dimension = D;
  const std::size_t rank;
  typedef BoundaryCondition3<dimension> condition_t;

  std::vector< boost::shared_ptr<condition_t> > conditions;

  typedef typename std::vector< boost::shared_ptr<condition_t> >::const_iterator const_iterator;

  const_iterator begin() const
  {
    return conditions.begin();
  }

  const_iterator end() const
  {
    return conditions.end();
  }

public:
  BoundaryConditionList(const std::size_t _rank) : rank(_rank)
  {
  }

  template<typename real_condition_t>
  void add(const real_condition_t& condition)
  {
    assert(condition.getRank() == rank);
    conditions.push_back(boost::shared_ptr<condition_t>(new real_condition_t(condition)));
  }

  virtual std::size_t getRank() const
  {
    return rank;
  }

  virtual bool applies(const int label) const
  {
    for(const_iterator condIter(begin()); condIter!=end(); ++condIter)
    {
      if ((*condIter)->applies(label))
        return true;
    }
    return false;
  }

  virtual int getPriority(const int label) const
  {
    std::size_t index = 0;

    for(const_iterator condIter(begin()); condIter!=end(); ++condIter)
    {
      --index;
      if ((*condIter)->applies(label))
        return index;
    }
    assert(false && "Tried to get priority of boundary condition with a label we don't have");
  }

  virtual Tensor<dimension> getValue(const vertex<dimension>& location, const int label) const
  {
    for(const_iterator condIter(begin()); condIter!=end(); ++condIter)
    {
      if ((*condIter)->applies(label))
        return (*condIter)->getValue(location, label);
    }

    assert(false && "Tried to get value for a boundary condition where it doesn't apply.");
    return Tensor<dimension>(rank);
  }
};

}

#endif
