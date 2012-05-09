#ifndef EXCAFE_SYMBOLIC_EXTRACTED_EXPRESSIONS_HPP
#define EXCAFE_SYMBOLIC_EXTRACTED_EXPRESSIONS_HPP

#include "symbolic_fwd.hpp"
#include "symbol.hpp"
#include "expr.hpp"
#include <utility>
#include <vector>
#include <set>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/multi_index/member.hpp>

namespace excafe
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

  std::size_t nextID;
  std::set<Symbol> representable;
  expr_map_t expressions;

  Symbol addExpression(const Expr& e, bool isRepresentable, bool topLevel);
  void getTopologicallySortedHelper(std::vector<Symbol>& result, std::set<Symbol>& seen, const Symbol& current) const;
  Symbol allocateSymbol(bool isRepresentable);

public:
  typedef expr_mapping_t value_type;
  typedef expr_map_t::index<symbol_tag>::type::iterator iterator;
  typedef expr_map_t::index<symbol_tag>::type::const_iterator const_iterator;

  ExtractedExpressions();

  const_iterator begin() const;
  const_iterator end() const;
  const_iterator find(const Symbol& s) const;

  Symbol addRepresentable(const Expr& e, bool topLevel = false);
  Symbol addUnrepresentable(const Expr& e, bool topLevel = false);
  void replaceExpression(const Symbol& s, const Expr& e);
  bool containsSymbol(const Symbol& s) const;
  bool isRepresentable(const Symbol& s) const;
  std::vector<Symbol> getTopologicallySorted() const;
};

}

}

#endif
