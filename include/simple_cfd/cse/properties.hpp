#ifndef SIMPLE_CFD_CSE_PROPERTIES_HPP
#define SIMPLE_CFD_CSE_PROPERTIES_HPP

namespace cfd
{

namespace cse
{

// Identities both the polynomial id and the term number
struct term_id       { typedef boost::edge_property_tag kind; };

// Polynomial id
struct polynomial_id { typedef boost::vertex_property_tag kind; };

// The value of the cube
struct term_cube     { typedef boost::vertex_property_tag kind; };

// The number of multiplies required to evaluate the cube/co-kernel
struct mul_count     { typedef boost::vertex_property_tag kind; };

// Is this a cube? Otherwise, it's a co-kernel.
struct is_cube       { typedef boost::vertex_property_tag kind; };

// Is this cube/co-kernel 1.0 or -1.0?
struct is_unit       { typedef boost::vertex_property_tag kind; };

// Is this cube a constant numeric value?
struct is_numeric    { typedef boost::vertex_property_tag kind; };

// Does this cube have a numeric coefficient?
struct has_coefficient    { typedef boost::vertex_property_tag kind; };


// This value represents a heuristic about which cubes to try to grow with first.
// It *must* have a unique value for each cube, hence the form (score, vertex_id).
struct cube_ordering { typedef boost::vertex_property_tag kind; };

}

}

#endif
