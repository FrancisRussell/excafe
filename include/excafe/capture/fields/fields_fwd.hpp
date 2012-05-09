#ifndef EXCAFE_CAPTURE_FIELDS_FIELDS_FWD_HPP
#define EXCAFE_CAPTURE_FIELDS_FIELDS_FWD_HPP

namespace excafe
{

class Scalar;
class Field;
class NamedField;
class FunctionSpace;
class Operator;
class TemporalIndex;
class BoundaryCondition;

namespace detail
{

// Function Space Related
class FunctionSpaceExpr;
class FunctionSpaceVisitor;
class FunctionSpaceMesh;
class FunctionSpaceUndefined;
class FunctionSpaceBinaryOperator;
class FunctionSpaceAddition;

class DiscreteExpr;
class DiscreteExprVisitor;

// Fields related
class DiscreteFieldExpr;
class DiscreteFieldUndefined;
class DiscreteFieldPersistent;
class DiscreteFieldZero;
class DiscreteFieldElementWise;
class DiscreteFieldTwoNorm;
class DiscreteFieldProjection;
class DiscreteFieldApplyBC;

// Scalar related
class ScalarExpr;
class ScalarLiteral;
class ScalarBinaryOperator;
class ScalarUndefined;

// Operator related
class OperatorExpr;
class OperatorAssembly;
class OperatorApplication;
class OperatorUndefined;
class OperatorAddition;
class OperatorApplyBC;

// Index related
class TemporalIndexExpr;
class TemporalIndexOffset;
class TemporalIndexValue;
class TemporalIndexSet;
class DiscreteIndexedScalar; 
class DiscreteIndexedField; 
class DiscreteIndexedOperator; 
template<typename discrete_object_tag> class IndexableValue;

// Solve related
class LinearSolve;

// Utility
class DiscreteExprContainerBuilder;
class DiscreteExprContainerBuilderFormVisitor;
}

}

#endif
