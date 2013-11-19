#include <excafe/codegen/codegen_fwd.hpp>
#include <excafe/codegen/product.hpp>
#include <excafe/util/lazy_copy.hpp>
#include <excafe/mp/float.hpp>
#include <excafe/numeric/cast.hpp>
#include <ostream>
#include <sstream>
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

std::string Product::negate(const std::string& exp, const bool sse) const
{
  if (sse)
  {
    assert(0 && "negate() should never be called for SSE");
  }
  else
  {
    return std::string("-") + exp;
  }
}

std::string Product::constructPositiveProduct(const exp_map_t& exps, const bool sse) const
{
  std::ostringstream stream;
  if (exps.empty())
  {
    if (sse)
      stream << "_mm_set_pd(1.0, 1.0)";
    else
      stream << "1.0";
  }
  else
  {
    if (sse)
    {
      // Count muliplicands
      int expCount = 0;
      BOOST_FOREACH(const exp_map_t::value_type& exp, exps)
      {
        expCount += exp.second;
      }

      // Generate calls to multiplies
      for(int i=1; i < expCount; ++i)
      {
        stream << "_mm_mul_pd(";
      }

      // Generate operands
      int opIndex = 0;
      BOOST_FOREACH(const exp_map_t::value_type& exp, exps)
      {
        for(int i=0; i<exp.second; ++i)
        {
          if (opIndex != 0)
            stream << ", ";

          stream << exp.first;

          if (opIndex != 0)
            stream << ")";

          ++opIndex;
        }
      }
    }
    else
    {
      bool firstDenominator = true;

      BOOST_FOREACH(const exp_map_t::value_type& exp, exps)
      {
        for(int i=0; i<exp.second; ++i)
        {
          if (!firstDenominator)
            stream << "*";
          else
            firstDenominator = false;

          stream << exp.first;
        }
      }
    }
  }

  return stream.str();
}

void Product::write(std::ostream& out, const bool sse) const
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

  if (!isUnitCoefficient || (sse && coefficient != 1.0))
  {
    std::ostringstream coeffStream;

    if (sse)
      coeffStream << "_mm_set_pd(" << coefficient << ", " << coefficient << ")";
    else
      coeffStream << coefficient;

    numerators.insert(std::make_pair(coeffStream.str(), 1));
  }

  std::string numeratorString = constructPositiveProduct(numerators, sse);
  std::string denominatorString = constructPositiveProduct(denominators, sse);

  std::ostringstream expStream;
  if (denominators.empty())
  {
    expStream << numeratorString;
  }
  else
  {
    if (sse)
      expStream << "_mm_div_pd(" << numeratorString << ", " << denominatorString << ")";
    else
      expStream << numeratorString << "/(" << denominatorString << ")";

  }

  if (coefficient < 0.0 && isUnitCoefficient && !sse)
    out << negate(expStream.str(), sse);
  else
    out << expStream.str();
}

void Product::write(std::ostream& out) const
{
  write(out, false);
}

void Product::writeSSE(std::ostream& out) const
{
  write(out, true);
}

}

}
