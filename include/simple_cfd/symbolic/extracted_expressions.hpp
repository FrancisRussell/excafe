#ifndef SIMPLE_CFD_SYMBOLIC_EXTRACTED_EXPRESSIONS_HPP
#define SIMPLE_CFD_SYMBOLIC_EXTRACTED_EXPRESSIONS_HPP

#include "symbolic_fwd.hpp"
#include "symbol.hpp"
#include "expr.hpp"
#include <utility>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/member.hpp>

namespace cfd
{

namespace symbolic
{

class ExtractedExpressions
{
private:
  struct symbol_tag {};
  struct expr_tag {};

  typedef std::pair<Symbol, Expr> expr_mapping_t;

  typedef boost::multi_index_container<
    expr_mapping_t,
    boost::multi_index::indexed_by<
      boost::multi_index::hashed_unique<
        boost::multi_index::tag<symbol_tag>,
        BOOST_MULTI_INDEX_MEMBER(expr_mapping_t, expr_mapping_t::first_type, first)>,
      boost::multi_index::hashed_non_unique<
        boost::multi_index::tag<expr_tag>,
        BOOST_MULTI_INDEX_MEMBER(expr_mapping_t, expr_mapping_t::second_type, second)>
    >
  > expr_map_t;

  expr_map_t representable;
  expr_map_t unrepresentable;

  static Symbol addExpression(expr_map_t& map, const Expr& e, bool newSymbol);

public:
  Symbol addRepresentable(const Expr& e, bool newSymbol = false);
  Symbol addUnrepresentable(const Expr& e, bool newSymbol = false);
};

}

}

#endif
