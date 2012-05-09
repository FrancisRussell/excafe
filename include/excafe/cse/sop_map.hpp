#ifndef EXCAFE_CSE_SOP_MAP_HPP
#define EXCAFE_CSE_SOP_MAP_HPP

#include "polynomial_index.hpp"
#include "sop.hpp"
#include <vector>
#include <map>
#include <cstddef>
#include <boost/utility.hpp>

namespace excafe
{

namespace cse
{

class SOPMap : public boost::noncopyable
{
private:
  std::size_t nextIndex;
  std::map<PolynomialIndex, SOP> map;

  const PolynomialIndex newIndex()
  {
    return PolynomialIndex(nextIndex++);
  }

public:
  typedef std::map<PolynomialIndex, SOP>::value_type value_type;
  typedef std::map<PolynomialIndex, SOP>::const_iterator const_iterator;
  typedef const_iterator iterator;

  SOPMap() : nextIndex(0)
  {
  }

  const_iterator begin() const;
  const_iterator end() const;

  PolynomialIndex reserveIndex();
  std::vector<PolynomialIndex> reserveIndices(const std::size_t count);
  PolynomialIndex addSOP(const SOP& sop);

  SOP& operator[](const PolynomialIndex& i);
  const SOP& operator[](const PolynomialIndex& i) const;
};

}

}

#endif
