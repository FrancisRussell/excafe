#ifndef SIMPLE_CFD_CSE_PROPERTIES_HPP
#define SIMPLE_CFD_CSE_PROPERTIES_HPP

namespace cfd
{

namespace cse
{

struct term_id       { typedef boost::edge_property_tag kind; };

struct term_cube     { typedef boost::vertex_property_tag kind; };
struct term_cokernel { typedef boost::vertex_property_tag kind; };
struct mul_count     { typedef boost::vertex_property_tag kind; };
struct is_cube       { typedef boost::vertex_property_tag kind; };

// This value represents a heuristic about which cubes to try to grow with first.
// It *must* have a unique value for each cube, hence the form (score, vertex_id).
struct cube_ordering { typedef boost::vertex_property_tag kind; };

}

}

#endif
