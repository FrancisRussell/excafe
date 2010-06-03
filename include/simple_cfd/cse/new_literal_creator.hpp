#ifndef SIMPLE_CFD_CSE_NEW_LITERAL_CREATOR_HPP
#define SIMPLE_CFD_CSE_NEW_LITERAL_CREATOR_HPP

#include "cse_fwd.hpp"

namespace cfd
{

namespace cse
{

class NewLiteralCreator
{
public:
  virtual unsigned getLiteralID(const PolynomialIndex& i) = 0;
};

}

}

#endif
