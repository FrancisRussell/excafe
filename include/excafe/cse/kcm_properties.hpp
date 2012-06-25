#ifndef EXCAFE_CSE_KCM_PROPERTIES_HPP
#define EXCAFE_CSE_KCM_PROPERTIES_HPP

#include "cube.hpp"
#include "polynomial_index.hpp"
#include <utility>

namespace excafe
{

namespace cse
{

struct VertexProperty
{
  // The value of the cube
  Cube term_cube;

  // Polynomial id
  PolynomialIndex polynomial_id;

  // This value represents a heuristic about which cubes to try to grow with first.
  // It *must* have a unique value for each cube, hence the form (score, vertex_id).
  typedef std::pair<int, unsigned> cube_ordering_t;
  cube_ordering_t cube_ordering;

  // The number of multiplies required to evaluate the cube/co-kernel
  int mul_count;

  // Is this a cube? Otherwise, it's a co-kernel.
  bool is_cube;

  // Is this cube/co-kernel 1.0 or -1.0?
  bool is_unit;

  // Is this cube a constant numeric value?
  bool is_numeric;

  // Does this cube have a numeric coefficient?
  bool has_coefficient;
};

struct EdgeProperty
{
  // The term number
  std::size_t term_id;
};

}

}

#endif
