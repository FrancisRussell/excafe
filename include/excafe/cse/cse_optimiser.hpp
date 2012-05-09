#ifndef EXCAFE_CSE_CSE_OPTIMISER_HPP
#define EXCAFE_CSE_CSE_OPTIMISER_HPP

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
#include <excafe/numeric/polynomial.hpp>
#include <excafe/numeric/cast.hpp>
#include <excafe/numeric/convert_expression.hpp>
#include <excafe/numeric/excafe_mapper.hpp>
#include <excafe/exception.hpp>
#include <excafe/symbolic/extracted_expressions.hpp>
#include "polynomial_index.hpp"
#include "factorised_expression_visitor.hpp"
#include "new_literal_creator.hpp"
#include "cube.hpp"
#include "sop.hpp"
#include "kcm.hpp"
#include "sop_builder.hpp"
#include "sop_map.hpp"

namespace excafe
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
  typedef symbolic::Symbol symbol_t;
  typedef typename NumericExpressionVisitor<symbolic::Expr>::float_t   float_t;
  typedef typename NumericExpressionVisitor<symbolic::Expr>::integer_t integer_t;
  typedef boost::variant<symbol_t, integer_t, float_t, polynomial_index_t> literal_t;

  SOPMap sops;
  symbolic::ExtractedExpressions extracted;
  boost::bimap<symbolic::Symbol, unsigned> originalIndices;
  boost::bimap<symbolic::Symbol, PolynomialIndex> sopSymbols;
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
      return iter->second;
    else
      CFD_EXCEPTION("Unknown literal ID.");
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

  class LiteralVisitor : public boost::static_visitor<symbolic::Expr>
  {
  private:
    const CSEOptimiser& parent;

  public:
    LiteralVisitor(const CSEOptimiser& _parent) : parent(_parent)
    {
    }

    result_type  operator()(const symbol_t& v) const
    {
      return v; 
    }
    
    result_type operator()(const integer_t& n) const
    {
      return n;
    }

    result_type operator()(const float_t& n) const
    {
      return n;
    }
    
    result_type operator()(const polynomial_index_t& index) const
    {
      const boost::bimap<symbolic::Symbol, PolynomialIndex>::right_const_iterator iter =
        parent.sopSymbols.right.find(index);

      if (iter != parent.sopSymbols.right.end())
        return iter->second;
      else
        CFD_EXCEPTION("Unable to map PolynomialIndex back to Symbol.");
    }
  };

  class FactorisedExpressionAdapter : public NumericExpressionVisitor<symbolic::Symbol>
  {
  private:
    const CSEOptimiser& parent;
    FactorisedExpressionVisitor<CSEOptimiser::variable_t>& child;
    std::map<symbolic::Symbol, PolynomialIndex> symbolIndices;

    PolynomialIndex getIndex(const symbolic::Symbol& s)
    {
      const std::map<symbolic::Symbol, PolynomialIndex>::const_iterator iter = symbolIndices.find(s);

      if (iter != symbolIndices.end())
      {
        return iter->second;
      }
      else
      {
        const PolynomialIndex newIndex(symbolIndices.size());
        symbolIndices.insert(std::make_pair(s, newIndex));
        return newIndex;
      }
    }

    bool getOriginalIndex(unsigned& index, const symbolic::Symbol& s)
    {
      const boost::bimap<symbolic::Symbol, unsigned>::left_const_iterator iter = 
        parent.originalIndices.left.find(s);
      
      if (iter != parent.originalIndices.left.end())
      {
        index = iter->second;
        return true;
      }
      else
      {
        return false;
      }
    }

    bool getPolynomialIndex(PolynomialIndex& index, const symbolic::Symbol& s)
    {
      const symbolic::ExtractedExpressions::const_iterator iter = 
        parent.extracted.find(s);

      if (iter != parent.extracted.end())
      {
        index = getIndex(s);
        return true;
      }
      else
      {
        return false;
      }
    }

  public:
    FactorisedExpressionAdapter(const CSEOptimiser& _parent, 
      FactorisedExpressionVisitor<CSEOptimiser::variable_t>& _child) :
      parent(_parent), child(_child)
    {
    }

    void visitExtracted()
    {
      using namespace excafe::symbolic;

      const ExtractedExpressions& extracted = parent.extracted;
      const std::vector<Symbol> sorted = extracted.getTopologicallySorted();

      BOOST_FOREACH(const Symbol& s, sorted)
      {
        const ExtractedExpressions::const_iterator iter = extracted.find(s);

        if (iter != extracted.end())
          iter->second.accept(*this);
        else
          CFD_EXCEPTION("Unable to find expression referenced in topological sort.");

        unsigned originalIndex;
        PolynomialIndex polynomialIndex;

        if (getOriginalIndex(originalIndex, s))
        {
          child.postOriginalTerm(originalIndex);
        }
        else if (getPolynomialIndex(polynomialIndex, s))
        {
          child.postFactorisedTerm(polynomialIndex);
        }
        else
        {
          CFD_EXCEPTION("Unable to determine what symbol maps back to.");
        }
      }
    }

    void visitConstant(const float_t& s)
    {
      child.visitConstant(s);
    }

    void visitConstant(const integer_t& s)
    {
      child.visitConstant(s);
    }

    void visitVariable(const variable_t& var)
    {
      using excafe::detail::ExcafeMapper;

      unsigned originalIndex;
      PolynomialIndex polynomialIndex;

      if (getOriginalIndex(originalIndex, var))
      {
        child.visitOriginalTerm(originalIndex);
      }
      else if (getPolynomialIndex(polynomialIndex, var))
      {
        child.visitFactorisedTerm(polynomialIndex);
      }
      else
      {
        const ExcafeMapper<CSEOptimiser::variable_t>& symbolMapper = ExcafeMapper<CSEOptimiser::variable_t>::instance();
        child.visitVariable(symbolMapper.getKey(var));
      }
    }

    void visitExponent(const int exponent)
    {
      child.visitExponent(exponent);
    }

    void visitAbsoluteValue()
    {
      child.visitAbsoluteValue();
    }

    void postSummation(const std::size_t nops)
    {
      child.postSummation(nops);
    }

    void postProduct(const std::size_t nops)
    {
      child.postProduct(nops);
    }
  };

  symbolic::Expr convertFromSOP(const SOP& sop)
  {
    symbolic::Expr result;

    for(SOP::const_iterator sopIter = sop.begin(); sopIter != sop.end(); ++sopIter)
    {
      symbolic::Expr term(1);

      for(Cube::const_iterator cubeIter = sopIter->begin(); cubeIter != sopIter->end(); ++cubeIter)
      {
        const LiteralVisitor literalVisitor(*this);
        const literal_t literal = getLiteral(cubeIter->first);
        const symbolic::Expr literalValue = boost::apply_visitor(literalVisitor, literal);

        term *= pow(literalValue, cubeIter->second);
      }
      result += term;
    }

    return result;
  }

public:
  template<typename InputIterator>
  CSEOptimiser(const InputIterator begin, const InputIterator end)
  {
    using excafe::detail::ExcafeMapper;
    using excafe::detail::convert_expression;
    using namespace symbolic;

    ExcafeMapper<variable_t>& symbolMapper = ExcafeMapper<variable_t>::instance();
    SOPBuilder<variable_t> builder(*this);
    std::size_t index = 0;

    for(InputIterator iter = begin; iter != end; ++iter)
    {
      const Expr original = convert_expression<Expr>(*iter, symbolMapper);
      const Expr rewritten = original.extractPolynomials(extracted);
      const Symbol symbol = extracted.addRepresentable(rewritten, true);

      originalIndices.left.insert(std::make_pair(symbol, index));
      ++index;
    }

    // Allocate indices for the extracted sub-expressions
    BOOST_FOREACH(const ExtractedExpressions::value_type& exprMapping, extracted)
    {
      const Symbol symbol = exprMapping.first;
      if (extracted.isRepresentable(symbol))
      {
        const PolynomialIndex polynomialIndex = sops.reserveIndex();
        sopSymbols.left.insert(std::make_pair(symbol, polynomialIndex));
        sops[polynomialIndex] = builder.getSOP(exprMapping.second);
      }
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

    // Create symbols for the newly constructed PolynomialIndex instances
    const Expr zero(0); 
    BOOST_FOREACH(const SOPMap::value_type& sopMapping, sops)
    {
      if (sopSymbols.right.find(sopMapping.first) == sopSymbols.right.end())
      {
        // We force top-level symbol creation to ensure symbol is
        // unique as we intend to redefine it later.
        const Symbol s = extracted.addRepresentable(zero, true);
        sopSymbols.left.insert(std::make_pair(s, sopMapping.first));
      }
    }

    // Copy all SOPs back into ExtractedExpressions
    BOOST_FOREACH(const SOPMap::value_type& sopMapping, sops)
    {
      const boost::bimap<Symbol, PolynomialIndex>::right_const_iterator sopSymIter = 
        sopSymbols.right.find(sopMapping.first);

      if (sopSymIter != sopSymbols.right.end())
      {
        const Symbol s = sopSymIter->second;
        const Expr e = convertFromSOP(sopMapping.second);
        extracted.replaceExpression(s, e);
      }
      else
      {
        CFD_EXCEPTION("Unable to map PolynomialIndex back to Symbol");
      }
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

  unsigned getLiteralID(const symbol_t& var)
  {
    return getLiteralIDTemplated(var);
  }

  unsigned getLiteralID(const integer_t& n)
  {
    return getLiteralIDTemplated(n);
  }

  unsigned getLiteralID(const float_t& n)
  {
    return getLiteralIDTemplated(n);
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
    const integer_t* const numeric = boost::get<integer_t>(&value);

    if (numeric == NULL)
    {
      return false;
    }
    else
    {
      const integer_t one(1);
      const integer_t minusOne(-1);
      return (*numeric) == one || (*numeric) == minusOne;
    }
  }

  bool isNumeric(const unsigned literal) const
  {
    const literal_t value = getLiteral(literal);
    return boost::get<integer_t>(&value) != NULL 
           || boost::get<float_t>(&value) != NULL;
  }

  bool isPolynomialIndex(const unsigned literal) const
  {
    const literal_t value = getLiteral(literal);
    return boost::get<polynomial_index_t>(&value) != NULL;
  }

  void accept(FactorisedExpressionVisitor<variable_t>& visitor) const
  {
    FactorisedExpressionAdapter adapter(*this, visitor);
    adapter.visitExtracted();
  }
};

}

}

#endif
