#ifndef SIMPLE_CFD_CSE_SOP_HPP
#define SIMPLE_CFD_CSE_SOP_HPP

#include <cstddef>
#include <map>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>
#include "cube.hpp"
#include <simple_cfd/exception.hpp>

namespace cfd
{

namespace cse
{

class SOP
{
public:
  typedef std::vector< std::pair<SOP, Cube> > kernel_set_t;

private:
  std::vector<unsigned> termNumbers;
  std::vector<Cube> cubes;

  void addCube(const Cube& cube)
  {
    const unsigned termNumber = cubes.size();
    addCube(termNumber, cube);
  }

  void addCube(const unsigned termNumber, const Cube& cube)
  {
    cubes.push_back(cube);
    termNumbers.push_back(termNumber);
  }

  void checkConsistent() const
  {
    if (termNumbers.size() != cubes.size())
      CFD_EXCEPTION("Inconsistency detected in SOP.");
  }

  static void addKernels(kernel_set_t& kernels, const unsigned i, const SOP& p, const Cube& d);

public:
  typedef std::vector<Cube>::value_type     value_type;
  typedef std::vector<Cube>::iterator       iterator;
  typedef std::vector<Cube>::const_iterator const_iterator;

  SOP()
  {
  }

  SOP(const Cube& cube)
  {
    addCube(cube);
  }

  SOP(const SOP& sop) : termNumbers(sop.termNumbers), cubes(sop.cubes)
  {
  }

  iterator begin()
  {
    return cubes.begin();
  }

  iterator end()
  {
    return cubes.end();
  }

  const_iterator begin() const
  {
    return cubes.begin();
  }

  const_iterator end() const
  {
    return cubes.end();
  }


  void append(const Cube& cube)
  {
    addCube(cube);
  }

  SOP operator/(const Cube& cube) const;

  Cube maxDivisor() const;
  std::map<unsigned, std::size_t> getLiteralUseCounts() const;
  kernel_set_t getKernels() const;

  template<typename literal_writer>
  void write(std::ostream& o, const literal_writer& writer) const
  {
    if (cubes.empty())
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
