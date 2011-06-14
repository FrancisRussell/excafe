#include <simple_cfd/symbolic/extracted_expressions.hpp>
#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>

namespace cfd
{

namespace symbolic
{

Symbol ExtractedExpressions::addExpression(expr_map_t& map, const Expr& e, const bool newSymbol)
{
  typedef expr_map_t::index<expr_tag>::type TaggedExprMap;
  const TaggedExprMap& taggedExprMap = map.get<expr_tag>();

  const TaggedExprMap::const_iterator exprIter = taggedExprMap.find(e);
  const bool notFound = (exprIter == taggedExprMap.end());

  if (newSymbol || notFound)
  {
    const Symbol newSym("extracted");
    const expr_mapping_t newEntry(newSym, (notFound ? e : exprIter->first));
    map.insert(newEntry);
    return newSym;
  }
  else
  {
    return exprIter->first;
  }
}

Symbol ExtractedExpressions::addRepresentable(const Expr& e, const bool newSymbol)
{
  return addExpression(representable, e, newSymbol);
}

Symbol ExtractedExpressions::addUnrepresentable(const Expr& e, const bool newSymbol)
{
  return addExpression(unrepresentable, e, newSymbol);
}

}

}
