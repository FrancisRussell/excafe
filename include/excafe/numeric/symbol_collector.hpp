#ifndef EXCAFE_NUMERIC_SYMBOL_COLLECTOR_HPP
#define EXCAFE_NUMERIC_SYMBOL_COLLECTOR_HPP

#include <set>
#include "expression_visitor.hpp"

namespace excafe
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
  void visitConstant(const typename parent_t::integer_t& value) {}
  void visitConstant(const typename parent_t::float_t& value) {}
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
