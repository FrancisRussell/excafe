#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FUNCTION_SPACE_VISITOR_HPP

#include "fields_fwd.hpp"

namespace cfd
{

namespace detail
{

class FunctionSpaceVisitor
{
public:
  virtual void visit(FunctionSpaceAddition& f) = 0;
  virtual void visit(FunctionSpaceMeshFunction& f) = 0;
  virtual void visit(FunctionSpaceUndefined& f) = 0;
  virtual ~FunctionSpaceVisitor() {}
};

}

}

#endif
