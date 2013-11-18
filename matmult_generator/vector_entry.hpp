#ifndef VECTOR_ENTRY_HPP
#define VECTOR_ENTRY_HPP

#include <iosfwd>

struct VectorEntry
{
  VectorEntry(const size_t _index) : index(_index)
  {
  }

  size_t index;

  bool operator<(const VectorEntry& e) const
  {
    return index < e.index;
  }
};

std::ostream& operator<<(std::ostream& o, const VectorEntry &e);

#endif
