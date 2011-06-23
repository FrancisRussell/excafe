#ifndef SIMPLE_CFD_NUMERIC_SYMBOL_COLLECTOR_HPP
#define SIMPLE_CFD_NUMERIC_SYMBOL_COLLECTOR_HPP

#include <set>
#include "expression_visitor.hpp"

namespace cfd
{

namespace detail
{

template<typename V>
class SymbolCollector : public NumericExpressionVisitor<V>
{
private:
  typedef V variable_t;
  typedef NumericExpressionVisitor<variable_t> parent_t;
  std::set<variable_t> symbols;

public:
  void visitConstant(const typename parent_t::value_t& value) {}
  void visitExponent(const int n) {}
  void visitAbsoluteValue() {}
  void postSummation(const std::size_t n) {}
  void postProduct(const std::size_t n) {}

  void visitVariable(const variable_t& v)
  {
    symbols.insert(v);
  }

  std::set<variable_t> getSymbols() const
  {
    return symbols;
  }
};

}

}

#endif
