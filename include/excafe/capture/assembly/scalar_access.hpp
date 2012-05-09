#ifndef EXCAFE_CAPTURE_ASSEMBLY_SCALAR_ACCESS_HPP
#define EXCAFE_CAPTURE_ASSEMBLY_SCALAR_ACCESS_HPP

#include <ostream>
#include <excafe/capture/fields/scalar.hpp>
#include <excafe/capture/fields/scalar_expr.hpp>

namespace excafe
{

namespace detail
{

class ScalarAccess
{
private:
  ScalarExpr::expr_ptr scalarExpr;

public:
  ScalarAccess(const Scalar& s) : scalarExpr(s.getExpr())
  {
  }

  bool operator==(const ScalarAccess& s) const
  {
    return scalarExpr == s.scalarExpr;
  }

  bool operator<(const ScalarAccess& s) const
  {
    return scalarExpr < s.scalarExpr;
  }

  ScalarExpr::expr_ptr getExpr() const
  {
    return scalarExpr;
  }

  void write(std::ostream& o) const
  {
    o << "scalar[" << scalarExpr << "]";
  }
};

std::ostream& operator<<(std::ostream& o, const ScalarAccess& c);

}

}

#endif
