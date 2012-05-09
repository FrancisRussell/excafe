#include <excafe/symbolic/extracted_expressions.hpp>
#include <excafe/symbolic/symbol.hpp>
#include <excafe/symbolic/expr.hpp>
#include <excafe/exception.hpp>
#include <sstream>
#include <vector>
#include <set>
#include <boost/foreach.hpp>

namespace excafe
{

namespace symbolic
{

ExtractedExpressions::ExtractedExpressions() : nextID(0)
{
}

ExtractedExpressions::const_iterator ExtractedExpressions::begin() const
{
  return expressions.get<symbol_tag>().begin();
}

ExtractedExpressions::const_iterator ExtractedExpressions::end() const
{
  return expressions.get<symbol_tag>().end();
}

ExtractedExpressions::const_iterator ExtractedExpressions::find(const Symbol& s) const
{
  return expressions.get<symbol_tag>().find(s);
}

Symbol ExtractedExpressions::allocateSymbol(const bool isRepresentable)
{
  std::ostringstream nameStream;
  nameStream << "extracted_" << nextID++;
  const Symbol symbol(nameStream.str());

  if (isRepresentable)
    representable.insert(symbol);

  return symbol;
}

Symbol ExtractedExpressions::addExpression(const Expr& e, const bool isRepresentable, const bool topLevel)
{
  typedef expr_map_t::index<expr_tag>::type TaggedExprMap;
  const TaggedExprMap& taggedExprMap = expressions.get<expr_tag>();

  const TaggedExprMap::const_iterator exprIter = taggedExprMap.find(e);
  const bool found = (exprIter != taggedExprMap.end());
  const bool createSymbol = !found || topLevel;

  if (createSymbol)
  {
    // Avoid creating new symbols for numeric values
    const Expr rhs = (is_exactly_a<Rational>(e) || !found) ? e : exprIter->first;

    const Symbol newSym = allocateSymbol(found || isRepresentable);
    const expr_mapping_t newEntry(newSym, rhs);
    expressions.insert(newEntry);
    return newSym;
  }
  else
  {
    return exprIter->first;
  }
}

Symbol ExtractedExpressions::addRepresentable(const Expr& e, const bool topLevel)
{
  return addExpression(e, true, topLevel);
}

Symbol ExtractedExpressions::addUnrepresentable(const Expr& e, const bool topLevel)
{
  return addExpression(e, false, topLevel);
}

void ExtractedExpressions::replaceExpression(const Symbol& s, const Expr& e)
{
  typedef expr_map_t::index<symbol_tag>::type TaggedExprMap;
  TaggedExprMap& taggedExprMap = expressions.get<symbol_tag>();
  const TaggedExprMap::iterator exprIter = taggedExprMap.find(s);

  if (exprIter != taggedExprMap.end())
  {
    const bool replaced = taggedExprMap.replace(exprIter, expr_mapping_t(s, e));
    assert(replaced);
  }
  else
  {
    CFD_EXCEPTION("Tried to replace symbol not present in ExtractedExpressions.");
  }
}

std::vector<Symbol> ExtractedExpressions::getTopologicallySorted() const
{
  std::vector<Symbol> result;
  std::set<Symbol> seen;

  BOOST_FOREACH(const value_type& mapping, *this)
    getTopologicallySortedHelper(result, seen, mapping.first);

  return result;
}

void ExtractedExpressions::getTopologicallySortedHelper(std::vector<Symbol>& result, 
                                                        std::set<Symbol>& seen, 
                                                        const Symbol& current) const
{
  if (seen.find(current) == seen.end())
    seen.insert(current);
  else
    return;
  
  const const_iterator iter = find(current);
  if (iter != end())
  {
    const std::set<Symbol> depends = iter->second.getSymbols();

    BOOST_FOREACH(const Symbol& symbol, depends)
      getTopologicallySortedHelper(result, seen, symbol);

    result.push_back(current);
  }
}

bool ExtractedExpressions::containsSymbol(const Symbol& s) const
{
  return find(s) != end();
}

bool ExtractedExpressions::isRepresentable(const Symbol& s) const
{
  assert(containsSymbol(s));
  return representable.find(s) != representable.end();
}

}

}
