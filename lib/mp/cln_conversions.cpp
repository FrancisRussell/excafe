#include <cln/cln.h>
#include <cln/float_io.h>
#include <sstream>
#include <excafe/mp/cln_conversions.hpp>
#include <excafe/mp/integer.hpp>
#include <excafe/mp/float.hpp>
#include <excafe/mp/rational.hpp>
#include <boost/lexical_cast.hpp>
#include <string>

namespace excafe
{

namespace mp
{

cln::cl_I  toCLN(const Integer& i)
{
  const std::string str = boost::lexical_cast<std::string>(i);
  return cln::cl_I(str.c_str());
}

cln::cl_RA toCLN(const Rational& r)
{
  const cln::cl_I numerator = toCLN(r.getNumerator());
  const cln::cl_I denominator = toCLN(r.getDenominator());
  return numerator/denominator;
}

cln::cl_F  toCLN(const Float& f)
{
  const std::string str = boost::lexical_cast<std::string>(f);
  return cln::cl_F(str.c_str());
}

Integer  fromCLN(const cln::cl_I& i)
{
  const std::string str = boost::lexical_cast<std::string>(i);
  return Integer(str.c_str());
}

Rational fromCLN(const cln::cl_RA& r)
{
  const Integer numerator = fromCLN(cln::numerator(r));
  const Integer denominator = fromCLN(cln::denominator(r));
  return Rational(numerator, denominator);
}

Float    fromCLN(const cln::cl_F& f)
{
  cln::cl_print_flags flags;
  flags.rational_base = 10;
  flags.default_float_format = cln::default_float_format;

  std::ostringstream stream;
  cln::print_float(stream, flags, f);
  const std::string str = stream.str();

  return Float(str.c_str());
}

}

}
