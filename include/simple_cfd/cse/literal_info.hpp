#ifndef SIMPLE_CFD_CSE_LITERAL_INFO_HPP
#define SIMPLE_CFD_CSE_LITERAL_INFO_HPP

#include <utility>

namespace cfd
{

namespace cse
{

class LiteralInfo
{
private:
  unsigned literal;
  bool reciprocal;

public:
  LiteralInfo(const unsigned _literal, const bool _reciprocal) : 
    literal(_literal), reciprocal(_reciprocal)
  {
  }

  LiteralInfo(const unsigned _literal) : 
    literal(_literal), reciprocal(false)
  {
  }

  unsigned getLiteral() const
  {
    return literal;
  }

  bool isReciprocal() const
  {
    return reciprocal;
  }

  bool operator<(const LiteralInfo& i) const
  {
    return std::make_pair(literal, reciprocal) 
           < std::make_pair(i.literal, i.reciprocal);
  }

  bool operator==(const LiteralInfo& i) const
  {
    return literal == i.literal 
           && reciprocal == i.reciprocal;
  }
};

}

}

#endif
