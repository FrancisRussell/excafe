#ifndef SIMPLE_CFD_CAPTURE_GENERIC_SYMBOL_HPP
#define SIMPLE_CFD_CAPTURE_GENERIC_SYMBOL_HPP

#include <cstddef>
#include <utility>
#include <ostream>
#include <string>

namespace cfd
{

namespace detail
{

class GenericSymbol
{
private:
  std::string name;

public:
  GenericSymbol(const std::string& _name) : name(_name)
  {
  }

  bool operator==(const GenericSymbol& s) const
  {
    return name == s.name;
  }

  bool operator<(const GenericSymbol& s) const
  {
    return name < s.name;
  }

  void write(std::ostream& o) const
  {
    o << name;
  }
};

std::ostream& operator<<(std::ostream& o, const GenericSymbol& c);

}

}

#endif
