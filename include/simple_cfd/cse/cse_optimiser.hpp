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
#include "expression_provider.hpp"
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

  class LiteralToC : public boost::static_visitor<std::string>
  {
  private:
    const CSEOptimiser& optimiser;
    ExpressionProvider<variable_t>& provider;
    bool declareSubterms;

  public:
    LiteralToC(const CSEOptimiser& _optimiser, 
               ExpressionProvider<variable_t>& _provider, 
               const bool _declareSubterms = false) : 
      optimiser(_optimiser), provider(_provider), declareSubterms(_declareSubterms)
    {
    }

    std::string operator()(const numeric_t& n) const
    {
      return boost::lexical_cast<std::string>(cfd::numeric_cast<double>(n));
    }

    std::string operator()(const variable_t& v) const
    {
      return provider.getRValue(v);
    }

    std::string operator()(const polynomial_index_t& i) const
    {
      if (optimiser.isOriginalExpression(i))
      {
        return provider.getLValue(optimiser.getOriginalIndex(i));
      }
      else
      {
        std::ostringstream result;
        result << (declareSubterms ? "const double " : "");
        result << "subterm_" << i.getIndex();
        return result.str();
      }
    }
  };

  std::string cubeToC(ExpressionProvider<variable_t>& provider, const Cube& cube) const
  {
    bool nullNumerator = true;
    bool nullDenominator = true;

    std::ostringstream numerator;
    std::ostringstream denominator;

    for(Cube::const_iterator cubeIter = cube.begin(); cubeIter != cube.end(); ++cubeIter)
    {
      bool& nullFlag = (cubeIter->second < 0 ? nullDenominator : nullNumerator);
      std::ostringstream& stream = (cubeIter->second < 0 ? denominator : numerator);

      const int absExp = std::abs(cubeIter->second);
      const literal_t literal = getLiteral(cubeIter->first);

      for(int i=0; i<absExp; ++i)
      {
        if (nullFlag)
          nullFlag = false;
        else
          stream << "*";
 
        stream << boost::apply_visitor(LiteralToC(*this, provider), literal);
      }
    }

    std::ostringstream result;

    if (nullNumerator)
      result << "1.0";
    else
      result << numerator.str();

    if (!nullDenominator)
      result << "/(" << denominator.str() << ")";

    return result.str();
  }

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

  void outputToC(ExpressionProvider<variable_t>& provider, std::ostream& out) const
  {
    const std::vector<PolynomialIndex> sorted = getTopologicallySorted();

    BOOST_FOREACH(const PolynomialIndex& index, sorted)
    {
      const SOP& sop = sops[index];

      out << LiteralToC(*this, provider, true)(index) << " = ";

      for(SOP::const_iterator sopIter = sop.begin(); sopIter != sop.end(); ++sopIter)
      {
        if (sopIter != sop.begin())
          out << " + ";
         
         const Cube& cube = *sopIter;
         out << cubeToC(provider, cube);
      }
      out << ";" << std::endl;
    }
  }
};

}

}

#endif
