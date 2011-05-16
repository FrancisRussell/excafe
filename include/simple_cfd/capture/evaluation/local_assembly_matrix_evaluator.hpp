#ifndef SIMPLE_CFD_CAPTURE_EVALUATION_LOCAL_ASSEMBLY_MATRIX_EVALUATOR_HPP
#define SIMPLE_CFD_CAPTURE_EVALUATION_LOCAL_ASSEMBLY_MATRIX_EVALUATOR_HPP

#include <memory>
#include <boost/shared_ptr.hpp>
#include "local_assembly_matrix_evaluator_impl.hpp"
#include <simple_cfd/capture/evaluation/expression_values.hpp>
#include <simple_cfd/local_assembly_matrix.hpp>

namespace cfd
{

namespace detail
{

template<std::size_t D>
class LocalAssemblyMatrixEvaluator
{
public:
  static const std::size_t dimension = D;

private:
  boost::shared_ptr< LocalAssemblyMatrixEvaluatorImpl<dimension> > impl;

public:
  LocalAssemblyMatrixEvaluator(std::auto_ptr< LocalAssemblyMatrixEvaluatorImpl<dimension> >& _impl) :
    impl(_impl.release())
  {
  }

  void evaluate(LocalAssemblyMatrix<dimension, double>& matrix,
                const std::size_t cid,
                const ExpressionValues<dimension>& values) const
  {
    impl->evaluate(matrix, cid, values);
  }
};

}

}

#endif
