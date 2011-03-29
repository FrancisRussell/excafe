#ifndef SIMPLE_CFD_CSE_CSE_OPTIMISER_HPP
#define SIMPLE_CFD_CSE_CSE_OPTIMISER_HPP

#include <map>
#include <vector>
#include <utility>
#include <cassert>
#include <boost/variant.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/exception.hpp>
#include "polynomial_index.hpp"
#include "new_literal_creator.hpp"
#include "cube.hpp"
#include "sop.hpp"
#include "kcm.hpp"
#include "sop_builder.hpp"
#include "sop_map.hpp"

namespace cfd
{

namespace cse
{

template<typename V>
class CSEOptimiser : public NewLiteralCreator
{
public:
  typedef V variable_t;

private:
  typedef PolynomialIndex polynomial_index_t;
  typedef typename NumericExpressionVisitor<variable_t>::value_t numeric_t;
  typedef boost::variant<variable_t, numeric_t, polynomial_index_t> literal_t;

  SOPMap sops;
  boost::bimap<literal_t, unsigned> literalNumbering;

  class TranslatedLiteralWriter
  {
  private:
    const CSEOptimiser& parent;

  public:
    TranslatedLiteralWriter(const CSEOptimiser& _parent) : parent(_parent)
    {
    }

    void operator()(std::ostream& o, const unsigned literal) const
    {
      o << parent.getLiteral(literal);
    }
  };

  template<typename L>
  unsigned getLiteralIDTemplated(const L& l)
  {
    const typename boost::bimap<literal_t, unsigned>::left_map::const_iterator iter = literalNumbering.left.find(l);

    if (iter != literalNumbering.left.end())
    {
      return iter->second;
    }
    else
    {
      const unsigned newID = literalNumbering.size();
      literalNumbering.left.insert(std::make_pair(l, newID));
      return newID;
    }
  }

  literal_t getLiteral(const unsigned id) const
  {
    const typename boost::bimap<literal_t, unsigned>::right_map::const_iterator iter = literalNumbering.right.find(id);

    if (iter != literalNumbering.right.end())
    {
      return iter->second;
    }
    else
    {
      CFD_EXCEPTION("Unknown literal ID.");
    }
  }

public:
  template<typename InputIterator>
  CSEOptimiser(const InputIterator begin, const InputIterator end)
  {
    std::vector<typename InputIterator::value_type> polynomials(begin, end);
    std::vector<PolynomialIndex> indices = sops.reserveIndices(polynomials.size());
    assert(polynomials.size() == indices.size());

    SOPBuilder<variable_t> builder(*this);
    for(std::size_t i=0; i<polynomials.size(); ++i)
      sops[indices[i]] = builder.getSOP(polynomials[i]);

    KCM kcm(*this);
    std::cout << "Cube count: " << kcm.numCubes() << std::endl;
    std::cout << "Co-kernel count: " << kcm.numCoKernels() << std::endl;
    std::cout << "Edge count: " << kcm.numEdges() << std::endl;
    
    while(kcm.factorise());

    const std::size_t numAdditions = kcm.numAdditions();
    const std::size_t numMultiplies = kcm.numMultiplies();

    std::cout << "Total floating point operations: " << numAdditions+numMultiplies;
    std::cout << " (" << numAdditions << " additions, " << numMultiplies << " multiplies)" << std::endl;

    std::cout << "Factorized SOPs:" << std::endl;
    BOOST_FOREACH(const SOPMap::value_type mapping, sops)
    {
      std::cout << mapping.first << ": ";
      write(std::cout, mapping.second);
      std::cout << std::endl;
    }
  }

  SOPMap& getSOPMap()
  {
    return sops;
  }

  const SOPMap& getSOPMap() const
  {
    return sops;
  }

  unsigned getLiteralID(const variable_t& var)
  {
    return getLiteralIDTemplated(var);
  }

  unsigned getLiteralID(const numeric_t& s)
  {
    return getLiteralIDTemplated(s);
  }

  unsigned getLiteralID(const polynomial_index_t& i)
  {
    return getLiteralIDTemplated(i);
  }

  void write(std::ostream& o, const Cube& c) const
  {
    c.write(o, TranslatedLiteralWriter(*this));
  }

  void write(std::ostream& o, const SOP& sop) const
  {
    sop.write(o, TranslatedLiteralWriter(*this));
  }
};

}

}

#endif
