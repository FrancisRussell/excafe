#ifndef SIMPLE_CFD_CSE_SOP_REWRITE_HPP
#define SIMPLE_CFD_CSE_SOP_REWRITE_HPP

#include <set>
#include <boost/foreach.hpp>
#include "cube.hpp"
#include "sop.hpp"

namespace cfd
{

namespace cse
{

class SOPRewrite
{
private:
  std::set<unsigned> removedTerms;
  std::set<Cube> addedCubes;

public:
  void addRemovedTerm(const unsigned literal)
  {
    removedTerms.insert(literal);
  }

  void addCube(const Cube& c)
  {
    addedCubes.insert(c);
  }

  SOP operator()(const SOP& sop) const
  {
    SOP result;

    for(SOP::const_iterator sopIter = sop.begin(); sopIter!=sop.end(); ++sopIter)
    {
      const std::size_t termNumber = sop.getTermNumber(sopIter);
      if (removedTerms.find(termNumber) == removedTerms.end())
      {
        result.addCube(termNumber, *sopIter);
      }
    }

    BOOST_FOREACH(const Cube& c, addedCubes)
    {
      result.append(c);
    }

    return result;
  }
};

}

}

#endif
