#include <ostream>
#include <cstddef>
#include <map>
#include <boost/foreach.hpp>
#include <simple_cfd/cse/sop.hpp>
#include <simple_cfd/cse/cube.hpp>
#include <simple_cfd/util/maybe.hpp>

namespace cfd
{

namespace cse
{

void SOP::addKernels(kernel_set_t& kernels, const unsigned i, const SOP& p, const Cube& d)
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

SOP SOP::operator/(const Cube& cube) const
{
  checkConsistent();

  SOP result;
  for(std::size_t i=0; i<cubes->size(); ++i)
  {
    const util::Maybe<Cube> dividedCube((*cubes)[i]/cube); 

    if (dividedCube.hasValue())
      result.addCube((*termNumbers)[i], dividedCube.value());
  }
  return result; 
}

Cube SOP::maxDivisor() const
{
  if (cubes->empty())
  {
    return Cube();
  }
  else
  {
    Cube result(*cubes->begin());
    BOOST_FOREACH(const Cube& c, *cubes)
    {
      result &= c;
    }
    return result;
  }
}

std::map<unsigned, std::size_t> SOP::getLiteralUseCounts() const
{
  std::map<unsigned, std::size_t> result;
  BOOST_FOREACH(const Cube& c, *cubes)
  {
    c.incrementUseCounts(result);
  }
  return result;
}

SOP::kernel_set_t SOP::getKernels() const
{
  kernel_set_t kernels;
  kernels.insert(kernels.end(), std::make_pair(*this, Cube()));
  addKernels(kernels, 0, *this, Cube());
  return kernels;
}

std::size_t SOP::getTermNumber(const const_iterator i) const
{
  const std::size_t offset = i - begin();
  return (*termNumbers)[offset];
}

std::ostream& operator<<(std::ostream& o, const SOP& sop)
{
  sop.write(o, detail::LiteralWriter());
  return o;
}

}

}
