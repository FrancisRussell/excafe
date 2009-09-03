#ifndef SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_FWD_HPP
#define SIMPLE_CFD_CAPTURE_FIELDS_FIELDS_FWD_HPP

namespace cfd
{

class Scalar;
class Field;
class NamedField;
class FunctionSpace;
class Operator;
class TemporalIndex;

namespace detail
{

// Function Space Related
class FunctionSpaceExpr;
class FunctionSpaceVisitor;
class FunctionSpaceMeshFunction;
class FunctionSpaceEmpty;
class FunctionSpaceBinaryOperator;
class FunctionSpaceAddition;

class DiscreteExprVisitor;

// Fields related
class DiscreteFieldExpr;
class DiscreteFieldUndefined;
class DiscreteFieldPersistent;
class DiscreteFieldZero;

// Scalar related
class ScalarExpr;
class ScalarLiteral;
class ScalarBinaryOperator;

// Operator related
class OperatorExpr;
class OperatorAssembly;
class OperatorApplication;
class OperatorUndefined;

// Index related
class TemporalIndexExpr;
class TemporalIndexOffset;
class TemporalIndexValue;
template<typename discrete_object_tag> class DiscreteObjectIndexed;
}

}

#endif
