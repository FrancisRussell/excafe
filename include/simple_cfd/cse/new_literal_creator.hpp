#ifndef SIMPLE_CFD_CSE_NEW_LITERAL_CREATOR_HPP
#define SIMPLE_CFD_CSE_NEW_LITERAL_CREATOR_HPP

#include <cstddef>
#include "cse_fwd.hpp"

namespace cfd
{

namespace cse
{

class NewLiteralCreator
{
public:
  virtual unsigned getLiteralID(const PolynomialIndex& i) = 0;
  virtual SOPMap& getSOPMap() = 0;
  virtual const SOPMap& getSOPMap() const = 0;
  virtual bool isUnit(const unsigned literal) const = 0;
  virtual bool isNumeric(const unsigned literal) const = 0;
};

}

}

#endif
