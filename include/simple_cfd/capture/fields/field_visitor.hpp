#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELD_VISITOR_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELD_VISITOR_HPP

#include "fields_fwd.hpp"

namespace cfd
{

namespace detail
{

class FieldVisitor
{
public:
  virtual void visit(FieldEmpty& e) = 0;
  virtual void visit(FieldPersistent& p) = 0;
  virtual ~FieldVisitor() {}
};

}

}


#endif
