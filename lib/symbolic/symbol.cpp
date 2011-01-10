#include <simple_cfd/symbolic/symbol.hpp>
#include <simple_cfd/symbolic/expr.hpp>
#include <simple_cfd/symbolic/number.hpp>
#include <boost/functional/hash.hpp>

namespace cfd
{

namespace symbolic
{

int Symbol::nextSerial = 0;

Symbol::Symbol(const std::string& _name) : name(_name), serial(nextSerial++)
{
}

std::size_t Symbol::nops() const
{
  return 0;
}

void Symbol::write(std::ostream& o) const
{
  o << name;
}

Expr Symbol::derivative(const Symbol& s) const
{
  if (serial == s.serial)
    return Expr(new Number(1));
  else
    return Expr(new Number(0));
}

bool Symbol::isNumber() const
{
  return false;
}

bool Symbol::operator==(const Symbol& s) const
{
  return serial == s.serial;
}

bool Symbol::operator<(const Symbol& s) const
{
  return serial < s.serial;
}

std::size_t Symbol::untypedHash() const
{
  return boost::hash<int>()(serial);
}

}

}
