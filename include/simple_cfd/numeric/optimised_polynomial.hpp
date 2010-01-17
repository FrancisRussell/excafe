#ifndef SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP
#define SIMPLE_CFD_NUMERIC_OPTIMISED_POLYNOMIAL_HPP

#include <simple_cfd_fwd.hpp>
#include <string>
#include <set>
#include <vector>
#include <utility>
#include <cstddef>

namespace cfd
{

class OptimisedPolynomial
{
private:
  typedef std::vector< std::pair<std::vector<std::size_t>, double> > coefficient_vec_t;
  const std::set<std::string> variables;
  coefficient_vec_t coefficients;
  mutable std::vector<double> paramData;

  std::vector<std::size_t> buildExponentVector(const Monomial<std::string>& m) const;
  double evaluate(const std::vector<double>& params) const;

public:
  OptimisedPolynomial(const Polynomial<std::string>& p);

  double operator()() const;
  double operator()(const double a) const;
  double operator()(const double a, const double b) const;
  double operator()(const double a, const double b, const double c) const;
};

}

#endif
