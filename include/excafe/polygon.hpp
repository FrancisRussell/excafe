#ifndef EXCAFE_POLYGON_HPP
#define EXCAFE_POLYGON_HPP

#include <cassert>
#include <cstddef>
#include <excafe_fwd.hpp>

namespace excafe
{

class Polygon
{
private:
  vertex<2> origin;
  std::size_t numSides;
  double radius;
  double rotation;

public:
  Polygon(const vertex<2>& _origin, const std::size_t _numSides, const double _radius, const double _rotation) :
    origin(_origin), numSides(_numSides), radius(_radius), rotation(_rotation)
  {
    assert(numSides > 2);
  }

  Polygon(const Polygon& p) : origin(p.origin), numSides(p.numSides), radius(p.radius), rotation(p.rotation)
  {
  }

  Polygon& operator=(const Polygon& p)
  {
    origin = p.origin;
    numSides = p.numSides;
    radius = p.radius;
    rotation = p.rotation;
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

  double getRotation() const
  {
    return rotation;
  }
};

}

#endif
