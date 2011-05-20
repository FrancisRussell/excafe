#include <ostream>
#include <cstddef>
#include <map>
#include <boost/foreach.hpp>
#include <simple_cfd/cse/sop.hpp>
#include <simple_cfd/cse/cube.hpp>
#include <simple_cfd/cse/literal_info.hpp>

namespace cfd
{

namespace cse
{

void SOP::addKernels(kernel_set_t& kernels, const unsigned i, const SOP& p, const Cube& d)
{
  typedef std::map<LiteralInfo, std::size_t> use_count_map;
  const use_count_map literalUseCounts = p.getLiteralUseCounts();

  BOOST_FOREACH(const use_count_map::value_type& useCountMapping, literalUseCounts)
  {
    const LiteralInfo jInfo = useCountMapping.first;
    const unsigned j = jInfo.getLiteral();
    if (j>=i && useCountMapping.second>1)
    {
      const Cube lj(j, jInfo.isReciprocal() ? -1 : 1);
      const SOP ft = p/lj;
      const Cube c = ft.maxDivisor();

      if (c.begin() == c.end() || c.begin()->first >= j)
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
  result.nextTermNumber = nextTermNumber;

  for(std::size_t i=0; i<cubes->size(); ++i)
  {
    const Cube& candidate = (*cubes)[i];

    if (candidate.contains(cube))
    {
      result.addCube((*termNumbers)[i], candidate-cube);
    }
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

std::map<LiteralInfo, std::size_t> SOP::getLiteralUseCounts() const
{
  std::map<LiteralInfo, std::size_t> result;
  BOOST_FOREACH(const Cube& c, *cubes)
  {
    c.incrementUseCounts(result);
  }
  return result;
}

std::size_t SOP::numAdditions(const NewLiteralCreator& creator) const
{
  std::size_t termCount = 0;
  bool hasNumeric = false;

  BOOST_FOREACH(const Cube& c, *cubes)
  {
    if (c.isNumeric(creator))
      hasNumeric = true;
    else
      ++termCount;
  }

  if (hasNumeric)
    ++termCount;

  return termCount > 0 ? termCount-1 : 0;
}

std::size_t SOP::numMultiplies(const NewLiteralCreator& creator) const
{
  std::size_t result = 0;
  BOOST_FOREACH(const Cube& c, *cubes)
  {
    result += c.numMultiplies(creator);
  }
  return result;
}

bool SOP::deleteTerm(const std::size_t termID)
{
  for(std::vector<unsigned>::iterator termNumIter = termNumbers->begin(); termNumIter!=termNumbers->end();
    ++termNumIter)
  {
    if (termID == *termNumIter)
    {
      const std::size_t offset = termNumIter - termNumbers->begin();
      termNumbers->erase(termNumbers->begin() + offset);
      cubes->erase(cubes->begin() + offset);
      return true;
    }
  }

  return false;
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
