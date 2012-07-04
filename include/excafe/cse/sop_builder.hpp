#ifndef EXCAFE_CSE_SOP_BUILDER_HPP
#define EXCAFE_CSE_SOP_BUILDER_HPP

#include "cse_fwd.hpp"
#include <excafe/numeric/factoriser.hpp>
#include <excafe/numeric/expression_visitor.hpp>
#include <excafe/numeric/cast.hpp>
#include <excafe/symbolic/expr.hpp>
#include <excafe/symbolic/symbol.hpp>
#include <excafe/exception.hpp>
#include <excafe/mp/integer.hpp>
#include <stack>
#include <set>

namespace excafe
{

namespace cse
{

template<typename V>
class SOPBuilder : public NumericExpressionVisitor<symbolic::Symbol>
{
private:
  typedef NumericExpressionVisitor<symbolic::Symbol> parent_t;
  typedef typename parent_t::variable_t variable_t;
  typedef typename parent_t::float_t    float_t;
  typedef typename parent_t::integer_t  integer_t;

  Factoriser factoriser;
  CSEOptimiser<V>& optimiser;

private:
  // SOPs don't enforce uniqueness so we use a set of cubes to represent an SOP
  // when we eliminate duplicate subexpressions.
  typedef std::multiset<Cube> sop_t;

  std::map<sop_t, unsigned> sopCache;
  std::stack<sop_t> stack;

  static SOP buildSOP(const sop_t& sopSet)
  {
    SOP result;
    BOOST_FOREACH(const Cube& c, sopSet)
      result.append(c);
    return result;
  }

  unsigned getCachedLiteral(const sop_t& sopSet)
  {
    const std::map<sop_t, unsigned>::const_iterator iter = sopCache.find(sopSet);

    if (iter != sopCache.end())
    {
      return iter->second;
    }
    else
    {
      SOP sop = buildSOP(sopSet);
      const PolynomialIndex index = optimiser.getSOPMap().addSOP(sop);
      const unsigned literal = optimiser.getLiteralID(index);
      sopCache.insert(std::make_pair(sopSet, literal));
      return literal;
    }
  }

  Cube popCube()
  {
    if (stack.empty())
      CFD_EXCEPTION("Tried to pop term from empty stack.");

    const sop_t terms = stack.top(); stack.pop();

    if (terms.empty())
      return Cube(optimiser.getLiteralID(integer_t(0)));
    else if (terms.size() == 1)
      return *terms.begin();
    else
      return Cube(getCachedLiteral(terms));
  }

  void pushCube(const Cube& c)
  {
    sop_t monomial;
    monomial.insert(c);
    stack.push(monomial);
  }

  void pushLiteral(const unsigned literalID)
  {
    pushCube(Cube(literalID));
  }

public:
  SOPBuilder(CSEOptimiser<V>& _optimiser) : optimiser(_optimiser)
  {
  }

  SOP getSOP(const symbolic::Expr& e)
  {
    assert(stack.empty());
    e.accept(*this);
    assert(stack.size() == 1);

    const SOP result = buildSOP(stack.top());
    stack.pop();
    return result;
  }

  void visitConstant(const float_t& s)
  {
    const unsigned literalID = optimiser.getLiteralID(s);
    pushLiteral(literalID);
  }

  void visitConstant(const integer_t& s)
  {
    typedef Factoriser::power_t power_t;

    const std::vector<power_t> factorMap = factoriser.factor(s);
    BOOST_FOREACH(const power_t& power, factorMap)
    {
      const unsigned literalID = optimiser.getLiteralID(power.first);
      pushLiteral(literalID);

      if (power.second != 1)
        visitExponent(power.second.toInt());
    }

    postProduct(factorMap.size());
  }

  void visitVariable(const symbolic::Symbol& var)
  {
    const unsigned literalID = optimiser.getLiteralID(var);
    pushLiteral(literalID);
  }

  void visitExponent(const int exponent)
  {
    Cube c = popCube();
    c *= exponent;
    pushCube(c);
  }

  void visitAbsoluteValue()
  {
    CFD_EXCEPTION("Cannot build a sum-of-products containing a modulus.");
  }

  void postSummation(const std::size_t nops)
  {
    if (stack.size() < nops)
      CFD_EXCEPTION("Too few operands on stack for summation");

    sop_t terms;
    for(std::size_t i=0; i<nops; ++i)
    {
      const sop_t subTerm = stack.top(); stack.pop();
      terms.insert(subTerm.begin(), subTerm.end());
    }
    stack.push(terms);
  }

  void postProduct(const std::size_t nops)
  {
    if (stack.size() < nops)
      CFD_EXCEPTION("Too few operands on stack for product");

    Cube result;

    for(std::size_t i=0; i<nops; ++i)
      result += popCube();

    pushCube(result);
  }
};

}

}

#endif
