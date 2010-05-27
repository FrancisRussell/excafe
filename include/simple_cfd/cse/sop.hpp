#ifndef SIMPLE_CFD_CSE_SOP_HPP
#define SIMPLE_CFD_CSE_SOP_HPP

#include <cstddef>
#include <map>
#include <vector>
#include <utility>
#include <boost/foreach.hpp>
#include "cube.hpp"
#include <simple_cfd/util/maybe.hpp>
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

  static void addKernels(kernel_set_t& kernels, const unsigned i, const SOP& p, const Cube& d)
  {
    typedef std::map<unsigned, std::size_t> use_count_map;
    const use_count_map literalUseCounts = p.getLiteralUseCounts();

    BOOST_FOREACH(const use_count_map::value_type& useCountMapping, literalUseCounts)
    {
      if (useCountMapping.second>1)
      {
        const unsigned j = useCountMapping.first;
        const Cube lj(j);
        const SOP ft = p/lj;
        const Cube c = ft.maxDivisor();

        if (c.isOne() || c.begin()->first >= j)
        {
          const SOP f1 = ft/c;
          const Cube d1 = merge(d, c, lj);

          kernels.insert(kernels.end(), std::make_pair(f1, d1));
          addKernels(kernels, j, f1, d1);
        }
      }
    }
  }

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

  SOP operator/(const Cube& cube) const
  {
    checkConsistent();

    SOP result;
    for(std::size_t i=0; i<cubes.size(); ++i)
    {
      const util::Maybe<Cube> dividedCube(cubes[i]/cube); 

      if (dividedCube.hasValue())
        result.addCube(termNumbers[i], dividedCube.value());
    }
    return result; 
  }

  Cube maxDivisor() const
  {
    if (cubes.empty())
    {
      return Cube();
    }
    else
    {
      Cube result(*cubes.begin());
      BOOST_FOREACH(const Cube& c, cubes)
      {
        result &= c;
      }
      return result;
    }
  }

  std::map<unsigned, std::size_t> getLiteralUseCounts() const
  {
    std::map<unsigned, std::size_t> result;
    BOOST_FOREACH(const Cube& c, cubes)
    {
      c.incrementUseCounts(result);
    }
    return result;
  }

  kernel_set_t getKernels() const
  {
    kernel_set_t kernels;
    kernels.insert(kernels.end(), std::make_pair(*this, Cube()));
    addKernels(kernels, 0, *this, Cube());
    return kernels;
  }

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
