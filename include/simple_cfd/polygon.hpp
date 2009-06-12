#ifndef SIMPLE_CFD_POLYGON_HPP
#define SIMPLE_CFD_POLYGON_HPP

#include <cassert>
#include <cstddef>
#include <simple_cfd_fwd.hpp>

namespace cfd
{

class Polygon
{
private:
  vertex<2> origin;
  std::size_t numSides;
  double radius;

public:
  Polygon(const vertex<2>& _origin, const std::size_t _numSides, const double _radius) :
    origin(_origin), numSides(_numSides), radius(_radius)
  {
    assert(numSides > 2);
  }

  Polygon(const Polygon& p) : origin(p.origin), numSides(p.numSides), radius(p.radius)
  {
  }

  Polygon& operator=(const Polygon& p)
  {
    origin = p.origin;
    numSides = p.numSides;
    radius = p.radius;
    return *this;
  }

  vertex<2> getOrigin() const
  {
    return origin;
  }

  std::size_t getNumSides() const
  {
    return numSides;
  }

  double getRadius() const
  {
    return radius;
  }
};

}

#endif
