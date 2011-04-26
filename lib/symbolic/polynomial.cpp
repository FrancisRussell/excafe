#include <simple_cfd/symbolic/polynomial.hpp>

namespace cfd
{

namespace symbolic
{


Polynomial::MonomialProduct Polynomial::Monomial::operator*(const Monomial& m) const
{
  return MonomialProduct(*this, m);
}

}

}
