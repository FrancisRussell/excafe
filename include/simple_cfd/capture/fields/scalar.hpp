#ifndef SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_SCALAR_HPP

#include "scalar_expr.hpp"
#include <boost/operators.hpp>

namespace cfd
{

class Scalar : boost::arithmetic<Scalar>
{
public:
  typedef detail::ScalarExpr::expr_ptr expr_ptr;

private:
  expr_ptr expr;

public:
  Scalar();
  Scalar(const double s);
  Scalar(detail::ScalarExpr* const expr);
  Scalar(detail::ScalarExpr::expr_ptr const expr);
  Scalar& operator+=(const Scalar& s);
  Scalar& operator-=(const Scalar& s);
  Scalar& operator*=(const Scalar& s);
  Scalar& operator/=(const Scalar& s);
  Scalar operator-() const;
  Scalar& operator=(const Scalar& s);
  Scalar operator<(const Scalar& s);
  Scalar operator<=(const Scalar& s);
  Scalar operator>(const Scalar& s);
  Scalar operator>=(const Scalar& s);
  Scalar operator==(const Scalar& s);
  expr_ptr getExpr() const;
};

}

#endif
