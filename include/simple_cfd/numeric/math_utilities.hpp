#ifndef SIMPLE_CFD_NUMERIC_MATH_UTILITIES_HPP
#define SIMPLE_CFD_NUMERIC_MATH_UTILITIES_HPP

#include <set>
#include <cstddef>
#include <string>
#include "numeric/numeric_fwd.hpp"

namespace cfd
{

class MathUtilities
{
private:
  static double jacobi_a_1_n(const double alpha, const double beta, const std::size_t n);
  static double jacobi_a_2_n(const double alpha, const double beta, const std::size_t n);
  static double jacobi_a_3_n(const double alpha, const double beta, const std::size_t n);
  static double jacobi_a_4_n(const double alpha, const double beta, const std::size_t n);

public:
  static Polynomial<std::string> jacobi(const double alpha, const double beta, const std::size_t n);
  static std::set<double> jacobi_roots(const double alpha, const double beta, const std::size_t n, const double epsilon = 1e-8);

  template<typename T>
  static T gcd(T u, T v)
  {
    u = (u < 0 ? -u : u);
    v = (v < 0 ? -v : v);
  
    if (u == 0 || v == 0)
      return u | v;
  
    unsigned shift;
    for(shift=0; ((u | v) & 1) == 0; ++shift)
    {
      u >>= 1;
      v >>= 1;
    }
  
    while ((u & 1) == 0)
      u >>= 1;
  
    do
    {
      while ((v & 1) == 0)
        v >>= 1;
  
      if (u > v)
      {
        const T tmp = u;
        u = v;
        v = tmp;
      }
  
      v -= u;
      v >>= 1;
    }
    while (v != 0);
  
    return u << shift;
  }

  template<typename T>
  static T lcm(const T a, const T b)
  {
    // Performing the division first reduces the size of the intermediate
    // value when gcd(a,b) > 1.

    return b*(a/gcd(a, b));
  }
};

}

#endif
