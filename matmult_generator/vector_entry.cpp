#include "vector_entry.hpp"

#include <ostream>
#include <sstream>

std::ostream& operator<<(std::ostream& o, const VectorEntry &e)
{
  std::ostringstream stream;
  stream << "v[" << e.index << "]";
  o << stream.str();
  return o;
}

