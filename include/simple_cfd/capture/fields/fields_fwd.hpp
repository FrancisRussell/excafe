#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_FWD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_FWD_HPP

namespace cfd
{

class Field;
class FunctionSpace;

namespace detail
{

// Function Space Related
class FunctionSpaceExpr;
class FunctionSpaceVisitor;
class FunctionSpaceMeshFunction;
class FunctionSpaceEmpty;
class FunctionSpaceBinaryOperator;
class FunctionSpaceAddition;


// Fields related
class FieldExpr;
class FieldVisitor;
class FieldEmpty;
class FieldPersistent;
}

}

#endif
