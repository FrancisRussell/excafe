#ifndef SIMPLE_CFD_CSE_SOP_HPP
#define SIMPLE_CFD_CSE_SOP_HPP

#include <cstddef>
#include <map>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>
#include "cse_fwd.hpp"
#include "cube.hpp"
#include "literal_info.hpp"
#include <simple_cfd/exception.hpp>
#include <simple_cfd/util/lazy_copy.hpp>

namespace cfd
{

namespace cse
{

class SOP
{
public:
  typedef std::vector<Cube> cokernel_set_t;

private:
  std::size_t nextTermNumber;
  util::LazyCopy< std::vector<unsigned> > termNumbers;
  util::LazyCopy< std::vector<Cube> > cubes;

  std::size_t addCube(const Cube& cube)
  {
    const std::size_t termNumber = nextTermNumber++;
    addCube(termNumber, cube);
    return termNumber;
  }

  void addCube(const unsigned termNumber, const Cube& cube)
  {
    cubes->push_back(cube);
    termNumbers->push_back(termNumber);
  }

  void checkConsistent() const
  {
    if (termNumbers->size() != cubes->size())
      CFD_EXCEPTION("Inconsistency detected in SOP.");
  }

  static void addCoKernels(cokernel_set_t& kernels, unsigned i, const SOP& p, const Cube& d);

public:
  typedef std::vector<Cube>::value_type     value_type;
  typedef std::vector<Cube>::iterator       iterator;
  typedef std::vector<Cube>::const_iterator const_iterator;

  SOP() : nextTermNumber(0)
  {
  }

  SOP(const Cube& cube) : nextTermNumber(0)
  {
    addCube(cube);
  }

  SOP(const SOP& sop) : nextTermNumber(sop.nextTermNumber), termNumbers(sop.termNumbers), cubes(sop.cubes)
  {
  }

  iterator begin()
  {
    return cubes->begin();
  }

  iterator end()
  {
    return cubes->end();
  }

  const_iterator begin() const
  {
    return cubes->begin();
  }

  const_iterator end() const
  {
    return cubes->end();
  }

  std::size_t size() const
  {
    return cubes->size();
  }

  std::size_t append(const Cube& cube)
  {
    return addCube(cube);
  }

  std::size_t numAdditions(const NewLiteralCreator& creator) const;
  std::size_t numMultiplies(const NewLiteralCreator& creator) const;
  bool deleteTerm(const std::size_t termID);

  SOP operator/(const Cube& cube) const;

  Cube maxDivisor() const;
  std::map<LiteralInfo, std::size_t> getLiteralUseCounts() const;
  cokernel_set_t getCoKernels() const;
  std::size_t getTermNumber(const const_iterator i) const;

  template<typename literal_writer>
  void write(std::ostream& o, const literal_writer& writer) const
  {
    if (cubes->empty())
      o << "0.0";

    for(SOP::const_iterator iter = begin(); iter != end(); ++iter)
    {
      if (iter != begin())
        o << " + ";

      iter->write(o, writer);
    }
  }
};

std::ostream& operator<<(std::ostream& o, const SOP& sop);

}

}

#endif
