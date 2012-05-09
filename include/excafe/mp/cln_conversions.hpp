#ifndef SIMPLE_CFD_MP_CLN_CONVERSIONS_HPP
#define SIMPLE_CFD_MP_CLN_CONVERSIONS_HPP

#include <cln/cln.h>
#include "mp_fwd.hpp"

namespace excafe
{

namespace mp
{

cln::cl_I  toCLN(const Integer& i);
cln::cl_RA toCLN(const Rational& r);
cln::cl_F  toCLN(const Float& f);

Integer  fromCLN(const cln::cl_I& i);
Rational fromCLN(const cln::cl_RA& r);
Float    fromCLN(const cln::cl_F& f);

}

}

#endif
