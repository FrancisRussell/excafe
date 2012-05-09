#ifndef EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_VISITOR_HPP
#define EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_VISITOR_HPP

#include "fields_fwd.hpp"

namespace excafe
{

namespace detail
{

class FunctionSpaceVisitor
{
public:
  virtual void visit(FunctionSpaceAddition& f) = 0;
  virtual void visit(FunctionSpaceMesh& f) = 0;
  virtual void visit(FunctionSpaceUndefined& f) = 0;
  virtual ~FunctionSpaceVisitor() {}
};

}

}

#endif
