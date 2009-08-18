#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_FWD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_FWD_HPP

namespace cfd
{

namespace detail
{

// Function Space Related
class FunctionSpace;
class FunctionSpaceExpr;
class FunctionSpaceVisitor;
class FunctionSpaceMeshFunction;
class FunctionSpaceEmpty;
class FunctionSpaceBinaryOperator;
class FunctionSpaceAddition;


// Fields related
class Field;
class FieldExpr;
class FieldVisitor;
class FieldEmpty;
class FieldPersistent;
}

}

#endif
