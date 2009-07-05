#ifndef SIMPLE_CFD_NUMERIC_QUADRATURE_HPP
#define SIMPLE_CFD_NUMERIC_QUADRATURE_HPP

#include <map>
#include <cstddef>

namespace cfd
{

class Quadrature
{
public:
  std::map<double, double> getGauss(const std::size_t n);
  std::map<double, double> getGaussRadau(const std::size_t n);
  std::map<double, double> getGaussLobatto(const std::size_t n);
};

}

#endif
