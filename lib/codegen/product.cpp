#include <excafe/codegen/codegen_fwd.hpp>
#include <excafe/codegen/product.hpp>
#include <excafe/util/lazy_copy.hpp>
#include <excafe/mp/float.hpp>
#include <excafe/numeric/cast.hpp>
#include <ostream>
#include <map>
#include <utility>
#include <boost/foreach.hpp>

namespace excafe
{

namespace codegen
{

Product Product::pow(const int exponent) const
{
  Product result;
  result.coefficient = mp::pow(coefficient, exponent);

  BOOST_FOREACH(const exp_map_t::value_type& exp, *exponents)
    (*result.exponents)[exp.first] = exp.second*exponent;

  return result;
}

Product& Product::operator*=(const Product& p)
{
  coefficient *= p.coefficient;
  BOOST_FOREACH(const exp_map_t::value_type& exp, *p.exponents)
    (*exponents)[exp.first] += exp.second;

  return *this;
}

void Product::write(std::ostream& out) const
{
  exp_map_t numerators, denominators;

  BOOST_FOREACH(const exp_map_t::value_type& exp, *exponents)
  {
    if (exp.second > 0)
      numerators.insert(exp);
    else
      denominators.insert(exp_map_t::value_type(exp.first, -exp.second));
  }

  const bool isUnitCoefficient = (coefficient == 1.0 || coefficient == -1.0);
  bool firstNumerator = true;

  if (numerators.empty() || !isUnitCoefficient)
  {
    firstNumerator = false;
    //out << std::setprecision(25) << coefficientAsDouble;
    out << coefficient;
  }
  else
  {
    out << (coefficient < 0.0 ? "-" : "");
  }

  BOOST_FOREACH(const exp_map_t::value_type& exp, numerators)
  {
    for(int i=0; i<exp.second; ++i)
    {
      if (!firstNumerator)
        out << "*";
      else
        firstNumerator = false;

      out << exp.first;
    }
  }

  if (!denominators.empty())
  {
    out << "/(";
    bool firstDenominator = true;

    BOOST_FOREACH(const exp_map_t::value_type& exp, denominators)
    {
      for(int i=0; i<exp.second; ++i)
      {
        if (!firstDenominator)
          out << "*";
        else
          firstDenominator = false;

        out << exp.first;
      }
    }

    out << ")";
  }
}

}

}
