#ifndef SIMPLE_CFD_CSE_CSE_OPTIMISER_HPP
#define SIMPLE_CFD_CSE_CSE_OPTIMISER_HPP

#include <map>
#include <vector>
#include <utility>
#include <cassert>
#include <ostream>
#include <string>
#include <sstream>
#include <boost/variant.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <boost/lexical_cast.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/numeric/cast.hpp>
#include <simple_cfd/exception.hpp>
#include "polynomial_index.hpp"
#include "factorised_expression_visitor.hpp"
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
  boost::bimap<unsigned, PolynomialIndex> originalIndices;
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

  bool isOriginalExpression(const PolynomialIndex& index) const
  {
    const boost::bimap<unsigned, PolynomialIndex>::right_map::const_iterator iter =
      originalIndices.right.find(index);

    return iter != originalIndices.right.end();
  }

  unsigned getOriginalIndex(const PolynomialIndex& index) const
  {
    const boost::bimap<unsigned, PolynomialIndex>::right_map::const_iterator iter =
      originalIndices.right.find(index);

    if (iter != originalIndices.right.end())
    {
      return iter->second;
    }
    else
    {
      CFD_EXCEPTION("Requested the original index for an expression that was created by the factoriser.");
    }
  }

  void getTopologicallySortedHelper(std::vector<PolynomialIndex>& result, 
                                    std::set<PolynomialIndex>& seen,
                                    const PolynomialIndex& current) const
  {
    if (seen.find(current) == seen.end())
      seen.insert(current);
    else
      return;

    const SOP& sop = sops[current];
    const std::map<LiteralInfo, std::size_t> useCounts = sop.getLiteralUseCounts();

    for(std::map<LiteralInfo, std::size_t>::const_iterator useIter = useCounts.begin(); 
        useIter != useCounts.end(); 
        ++useIter)
    {
      const unsigned literal = useIter->first.getLiteral();
      if (isPolynomialIndex(literal))
      {
        const PolynomialIndex index = boost::get<polynomial_index_t>(getLiteral(literal));
        getTopologicallySortedHelper(result, seen, index);
      }
    }

    result.push_back(current);
  }

  std::vector<PolynomialIndex> getTopologicallySorted() const
  {
    std::vector<PolynomialIndex> result;
    std::set<PolynomialIndex> seen;

    typedef boost::bimap<unsigned, PolynomialIndex>::left_map::value_type index_mapping_t;
    BOOST_FOREACH(const index_mapping_t& mapping, originalIndices.left)
    {
      getTopologicallySortedHelper(result, seen, mapping.second);
    }

    return result;
  }

  class LiteralVisitor : public boost::static_visitor<void>
  {
  private:
    const CSEOptimiser& parent;
    FactorisedExpressionVisitor<variable_t>& visitor;

  public:
    LiteralVisitor(const CSEOptimiser& _parent, FactorisedExpressionVisitor<variable_t>& _visitor) :
      parent(_parent), visitor(_visitor)
    {
    }

    void operator()(const variable_t& v) const
    {
      visitor.visitVariable(v);
    }
    
    void operator()(const numeric_t& n) const
    {
      visitor.visitConstant(n);
    }
    
    void operator()(const polynomial_index_t& index) const
    {
      if (parent.isOriginalExpression(index))
        visitor.visitOriginalTerm(parent.getOriginalIndex(index));
      else
        visitor.visitFactorisedTerm(index);
    }
  };

public:
  template<typename InputIterator>
  CSEOptimiser(const InputIterator begin, const InputIterator end)
  {
    const std::vector<typename InputIterator::value_type> polynomials(begin, end);
    const std::vector<PolynomialIndex> indices = sops.reserveIndices(polynomials.size());
    assert(polynomials.size() == indices.size());

    SOPBuilder<variable_t> builder(*this);
    for(std::size_t i=0; i<polynomials.size(); ++i)
    {
      originalIndices.left.insert(std::make_pair(i, indices[i]));
      sops[indices[i]] = builder.getSOP(polynomials[i]);
    }

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

  bool isUnit(const unsigned literal) const
  {
    const literal_t value = getLiteral(literal);
    const numeric_t* const numeric = boost::get<numeric_t>(&value);

    if (numeric == NULL)
    {
      return false;
    }
    else
    {
      static const numeric_t one(1);
      static const numeric_t minusOne(-1);
      return (*numeric) == one || (*numeric) == minusOne;
    }
  }

  bool isNumeric(const unsigned literal) const
  {
    const literal_t value = getLiteral(literal);
    return (boost::get<numeric_t>(&value) != NULL);
  }

  bool isPolynomialIndex(const unsigned literal) const
  {
    const literal_t value = getLiteral(literal);
    return (boost::get<polynomial_index_t>(&value) != NULL);
  }

  void accept(FactorisedExpressionVisitor<variable_t>& visitor) const
  {
    const std::vector<PolynomialIndex> sorted = getTopologicallySorted();

    BOOST_FOREACH(const PolynomialIndex& index, sorted)
    {
      const SOP& sop = sops[index];
      for(SOP::const_iterator sopIter = sop.begin(); sopIter != sop.end(); ++sopIter)
      {
        for(Cube::const_iterator cubeIter = sopIter->begin(); cubeIter != sopIter->end(); ++cubeIter)
        {
          const LiteralVisitor literalVisitor(*this, visitor);
          const literal_t literal = getLiteral(cubeIter->first);
          boost::apply_visitor(literalVisitor, literal);

          if (cubeIter->second != 1)
            visitor.visitExponent(cubeIter->second);
        }
        visitor.postProduct(sopIter->size());
      }
      visitor.postSummation(sop.size());

      if (isOriginalExpression(index))
        visitor.postOriginalTerm(getOriginalIndex(index));
      else
        visitor.postFactorisedTerm(index);
    }
  }
};

}

}

#endif
