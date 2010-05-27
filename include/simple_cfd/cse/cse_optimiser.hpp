#ifndef SIMPLE_CFD_CSE_CSE_OPTIMISER_HPP
#define SIMPLE_CFD_CSE_CSE_OPTIMISER_HPP

#include <map>
#include <vector>
#include <utility>
#include <boost/variant.hpp>
#include <boost/foreach.hpp>
#include <boost/bimap.hpp>
#include <simple_cfd/numeric/polynomial.hpp>
#include <simple_cfd/exception.hpp>
#include "cube.hpp"
#include "sop.hpp"

namespace cfd
{

namespace cse
{

template<typename V>
class CSEOptimiser
{
public:
  typedef V variable_t;

private:
  typedef Polynomial<variable_t> polynomial_t;
  typedef typename Polynomial<variable_t>::monomial_t monomial_t;
  typedef typename polynomial_t::internal_value_t numeric_t;
  typedef boost::variant<variable_t, numeric_t> literal_t;
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

  Cube buildCube(const numeric_t& coefficient, const monomial_t& monomial)
  {
    std::vector< std::pair<unsigned, unsigned> > exponentSeq;
    exponentSeq.reserve(monomial.size()+1);

    if (coefficient != 1.0)
      exponentSeq.push_back(std::make_pair(getLiteralID(coefficient), 1u));

    BOOST_FOREACH(const typename monomial_t::value_type& exponent, monomial)
    {
      exponentSeq.push_back(std::make_pair(getLiteralID(exponent.first), exponent.second));
    }

    return Cube(exponentSeq.begin(), exponentSeq.end());
  }

  SOP buildSOP(const polynomial_t& p)
  {
    typedef std::pair<monomial_t, numeric_t> monomial_coefficient;

    SOP sop;
    BOOST_FOREACH(const monomial_coefficient& monomialCoefficient, p)
    {
      sop.append(buildCube(monomialCoefficient.second, monomialCoefficient.first));
    }
    return sop;
  }

public:
  template<typename InputIterator>
  CSEOptimiser(const InputIterator begin, const InputIterator end)
  {
    std::vector<polynomial_t> polynomials(begin, end);
    std::vector<SOP> sops;

    BOOST_FOREACH(const polynomial_t& p, polynomials)
    {
      sops.push_back(buildSOP(p));
    }

    BOOST_FOREACH(const SOP& sop, sops)
    {
      const SOP::kernel_set_t kernels = sop.getKernels();
      
      BOOST_FOREACH(const SOP::kernel_set_t::value_type kernel, kernels)
      {
        std::cout << "kernel: ";
        write(std::cout, kernel.first);
        std::cout << ", co-kernel: ";
        write(std::cout, kernel.second);
        std::cout << std::endl;
      }
    }
  }

  unsigned getLiteralID(const variable_t& var)
  {
    return getLiteralIDTemplated(var);
  }

  unsigned getLiteralID(const numeric_t& s)
  {
    return getLiteralIDTemplated(s);
  }

  void write(std::ostream& o, const Cube& c) const
  {
    c.write(o, TranslatedLiteralWriter(*this));
  }

  void write(std::ostream& o, const SOP& sop) const
  {
    sop.write(o, TranslatedLiteralWriter(*this));
  }
};

}

}

#endif
