#ifndef SIMPLE_CFD_SYMBOLIC_EXPR_HPP
#define SIMPLE_CFD_SYMBOLIC_EXPR_HPP

#include <boost/shared_ptr.hpp>
#include <map>
#include "symbolic_fwd.hpp"

namespace cfd
{

namespace symbolic
{

class Expr
{
public:
  typedef boost::shared_ptr<const Basic> ref_t;

private:
  ref_t expr;

public:
  typedef std::map<Symbol, Expr> subst_map;

  explicit Expr(Basic* e);
  explicit Expr(ref_t& e);
  Expr();
  Expr(const double s);
  Expr& operator=(const Expr& e);
  bool operator<(const Expr& e) const;
  bool operator==(const Expr& e) const;
  bool operator!=(const Expr& e) const;
  void write(std::ostream& o) const;
  std::size_t hashValue() const;
  Expr derivative(const Symbol& s) const;
  Expr simplify() const;
  Expr integrate(const Symbol& s) const;
  const Basic& internal() const;
  Expr subs(const subst_map& map) const;
  virtual void accept(Visitor& v) const;
  virtual Expr expand() const;
};

std::size_t hash_value(const Expr& e);
Expr operator+(const Expr& a, const Expr& b);
Expr operator-(const Expr& a, const Expr& b);
Expr operator*(const Expr& a, const Expr& b);
Expr operator/(const Expr& a, const Expr& b);
Expr pow(const Expr& e, int power);
std::ostream& operator<<(std::ostream& o, const Expr& e);

}

}


#endif
