#ifndef EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_EXPR_HPP 
#define EXCAFE_CAPTURE_FIELDS_FUNCTION_SPACE_EXPR_HPP 

#include <boost/shared_ptr.hpp>
#include "fields_fwd.hpp"

namespace excafe
{

namespace detail
{

class FunctionSpaceExpr
{
public:
  typedef boost::shared_ptr<FunctionSpaceExpr> expr_ptr;

  virtual void accept(FunctionSpaceVisitor& visitor) = 0;
  virtual ~FunctionSpaceExpr() {}
};

}

}

#endif
