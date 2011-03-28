#ifndef SIMPLE_CFD_CSE_SOP_BUILDER_HPP
#define SIMPLE_CFD_CSE_SOP_BUILDER_HPP

#include "cse_fwd.hpp"
#include <simple_cfd/numeric/factoriser.hpp>
#include <simple_cfd/numeric/expression_visitor.hpp>
#include <simple_cfd/exception.hpp>
#include <stack>
#include <set>

#include <cln/io.h>
#include <cln/ring.h>
#include <cln/integer_ring.h>

namespace cfd
{

namespace cse
{

template<typename V>
class SOPBuilder : public NumericExpressionVisitor<V>
{
private:
  typedef NumericExpressionVisitor<V> parent_t;
  typedef typename parent_t::variable_t variable_t;
  typedef typename parent_t::value_t    value_t;

  Factoriser factoriser;
  CSEOptimiser<variable_t>& optimiser;

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
    {
      return Cube(optimiser.getLiteralID(value_t(0)));
    }
    else if (terms.size() == 1)
    {
      return *terms.begin();
    }
    else
    {
      return Cube(getCachedLiteral(terms));
    }
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
  SOPBuilder(CSEOptimiser<variable_t>& _optimiser) : optimiser(_optimiser)
  {
  }

  void visitConstant(const value_t& s)
  {
    // Checks if s is an integer.
    if (cln::instanceof(s, cln::cl_I_ring))
    {
      typedef Factoriser::power_t power_t;

      const std::vector<power_t> factorMap = factoriser.factor(cln::the<cln::cl_I>(s));
      BOOST_FOREACH(const power_t& power, factorMap)
      {
        const unsigned literalID = optimiser.getLiteralID(power.first);
        pushLiteral(literalID);

        if (power.second != 1)
          visitExponent(cln::cl_I_to_int(power.second));
      }

      postProduct(factorMap.size());
    }
    else
    {
      const unsigned literalID = optimiser.getLiteralID(s);
      pushLiteral(literalID);
    }
  }

  void visitVariable(const variable_t& var)
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

  SOP getResult() const
  {
    assert(stack.size() == 1);
    return buildSOP(stack.top());
  }
};

}

}

#endif
