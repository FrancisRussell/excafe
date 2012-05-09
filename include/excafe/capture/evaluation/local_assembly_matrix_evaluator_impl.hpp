#ifndef EXCAFE_CAPTURE_EVALUATION_LOCAL_ASSEMBLY_MATRIX_EVALUATOR_IMPL_HPP
#define EXCAFE_CAPTURE_EVALUATION_LOCAL_ASSEMBLY_MATRIX_EVALUATOR_IMPL_HPP

namespace excafe
{

namespace detail
{

template<std::size_t D>
class LocalAssemblyMatrixEvaluatorImpl
{
public:
  static const std::size_t dimension = D;
  typedef ScalarPlaceholder::expression_t expression_t;

  virtual void evaluate(LocalAssemblyMatrix<dimension, double>& matrix,
                        std::size_t cid,
                        const ExpressionValues<dimension>& values) const = 0;

  virtual ~LocalAssemblyMatrixEvaluatorImpl() {}
};

}

}

#endif
