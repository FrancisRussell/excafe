#ifndef SIMPLE_CFD_FE_MATRIX_HPP
#define SIMPLE_CFD_FE_MATRIX_HPP

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <boost/foreach.hpp>
#include "vertex.hpp"
#include "dof_map.hpp"
#include "discrete_field.hpp"
#include "mesh.hpp"
#include "numeric/matrix.hpp"
#include "numeric/functional.hpp"
#include "numeric/sparsity_pattern.hpp"
#include "quadrature_points.hpp"
#include "local_assembly_matrix.hpp"
#include "exception.hpp"
#include "capture/forms/bilinear_form_integral_sum.hpp"
#include "capture/forms/basis_finder.hpp"
#include "capture/assembly/assembly_helper.hpp"
#include "capture/assembly/scalar_placeholder.hpp"
#include "capture/assembly/scalar_placeholder_evaluator.hpp"
#include "capture/evaluation/local_assembly_matrix_evaluator.hpp"

namespace cfd
{

template<std::size_t D>
class DiscreteOperator
{
private:
  static const std::size_t dimension = D;

  typedef vertex<dimension>                 vertex_type;
  typedef FiniteElement<dimension>          finite_element_t;
  typedef typename DofMap<dimension>::dof_t dof_t;

  // Useful typedefs for expression capture
  typedef detail::ScalarPlaceholder::expression_t              expression_t;
  typedef detail::ScalarPlaceholder::optimised_expression_t    optimised_expression_t;
  typedef detail::LocalAssemblyMatrix<dimension, expression_t> local_matrix_t;
  typedef detail::LocalAssemblyMatrix<dimension, double>       evaluated_local_matrix_t;
  typedef detail::LocalAssemblyMatrixEvaluator<dimension>      cell_integral_t;

  const DofMap<dimension> rowMappings;
  const DofMap<dimension> colMappings;
  PETScMatrix matrix;

  static SparsityPattern createSparsityPattern(const DofMap<dimension>& rowMappings,
                                               const DofMap<dimension>& colMappings)
  {
    const unsigned rowDofs = rowMappings.getDegreesOfFreedomCount();
    const unsigned colDofs = colMappings.getDegreesOfFreedomCount();

    SparsityPattern pattern(rowDofs, colDofs);
    const Mesh<dimension>& m(rowMappings.getMesh());

    const std::set<const finite_element_t*> rowElements(rowMappings.getFiniteElements());
    const std::set<const finite_element_t*> colElements(colMappings.getFiniteElements());

    for(typename Mesh<dimension>::global_iterator cellIter(m.global_begin(dimension)); cellIter!=m.global_end(dimension); ++cellIter)
    {
      for(typename std::set<const finite_element_t*>::const_iterator rowElemIter(rowElements.begin()); rowElemIter!=rowElements.end(); ++rowElemIter)
      {
        for(typename std::set<const finite_element_t*>::const_iterator colElemIter(colElements.begin()); colElemIter!=colElements.end(); ++colElemIter)
        {
          for(unsigned rowDof=0; rowDof < (*rowElemIter)->spaceDimension(); ++rowDof)
          {
            for(unsigned colDof=0; colDof < (*colElemIter)->spaceDimension(); ++colDof)
            {
              const int rowIndex =
                rowMappings.getGlobalIndexWithMissingAsNegative(dof_t(*rowElemIter,
                cellIter->getIndex(), rowDof)); 
              const int colIndex =
                colMappings.getGlobalIndexWithMissingAsNegative(dof_t(*colElemIter,
                cellIter->getIndex(), colDof));

              if (rowIndex >= 0 && colIndex >= 0)
                pattern.insert(rowIndex, colIndex); 
            }
          }
        }
      }
    }

    return pattern;
  }

  template<typename value_map_t>
  value_map_t evaluatePlaceholders(const Scenario<dimension>& scenario,
    const detail::ExpressionValues<dimension>& expressionValues, const std::size_t cid, 
    const std::set<detail::ScalarPlaceholder>& placeholders) const
  {
    using namespace detail;

    const ScalarPlaceholderEvaluator<dimension> evaluator(scenario, expressionValues, cid);
    value_map_t values;

    BOOST_FOREACH(const ScalarPlaceholder& placeholder, placeholders)
    {
      values.bind(placeholder, evaluator(placeholder));
    }

    return values;
  }

public:
  DiscreteOperator(const DofMap<dimension>& _rowMappings, const DofMap<dimension>& _colMappings) :
          rowMappings(_rowMappings), colMappings(_colMappings), matrix(createSparsityPattern(rowMappings, colMappings))
  {
    // Make sure the degree-of-freedom maps are defined on the same mesh. Calculating a
    // sparsity pattern doesn't make sense otherwise.
    assert(&rowMappings.getMesh() == &colMappings.getMesh());
  }

  DiscreteOperator(const DofMap<dimension>& _rowMappings, const DofMap<dimension>& _colMappings, const PETScMatrix& m) :
          rowMappings(_rowMappings), colMappings(_colMappings), matrix(m)
  {
    assert(&rowMappings.getMesh() == &colMappings.getMesh());
  }

  DiscreteOperator(const DiscreteOperator& m) : rowMappings(m.rowMappings), colMappings(m.colMappings), matrix(m.matrix)
  {
  }

  DofMap<dimension> getRowMappings() const
  {
    return rowMappings;
  }

  DofMap<dimension> getColMappings() const
  {
    return colMappings;
  }

  DiscreteOperator& operator=(const DiscreteOperator& f)
  {
    assert(rowMappings == f.rowMappings);
    assert(colMappings == f.colMappings);
    matrix = f.matrix;
    return *this;
  }

  void addValues(const std::size_t cid, const detail::LocalAssemblyMatrix<dimension, double>& localMatrix)
  {
    //std::cout << "Adding local assembly matrix for cell " << cid << ":" << localMatrix << std::endl;
    const std::vector<dof_t> testDofs(localMatrix.getTestDofs(cid));
    const std::vector<dof_t> trialDofs(localMatrix.getTrialDofs(cid));
    addValues(testDofs.size(), trialDofs.size(), &testDofs[0], &trialDofs[0], localMatrix.data());
  }

  void addValues(const unsigned rows, const unsigned cols, const dof_t* rowDofs, const dof_t* colDofs, const double* block)
  {
    std::vector<int> rowIndices(rows);
    std::vector<int> colIndices(cols);

    for(unsigned row=0; row<rows; ++row)
      rowIndices[row] = rowMappings.getGlobalIndex(rowDofs[row]);

    for(unsigned col=0; col<cols; ++col)
      colIndices[col] = colMappings.getGlobalIndex(colDofs[col]);

    matrix.addValues(rows, cols, &rowIndices[0], &colIndices[0], block);
  }

  void assembleFromOptimisedLocalMatrix(const Scenario<dimension>& scenario, 
    const detail::ExpressionValues<dimension>& values,
    const std::map<MeshEntity, cell_integral_t> localMatrices,
    const MeshFunction<bool>& subDomain)
  {
    using namespace detail;

    const std::set<const finite_element_t*> trialElements(colMappings.getFiniteElements());
    const std::set<const finite_element_t*> testElements(rowMappings.getFiniteElements());
    const Mesh<dimension>& m = scenario.getMesh();
    const std::size_t entityDimension = subDomain.getDimension();
    evaluated_local_matrix_t localMatrix(testElements, trialElements);

    std::cout << "Starting assembly...." << std::flush;
    for(typename Mesh<dimension>::global_iterator eIter(m.global_begin(entityDimension)); eIter != m.global_end(entityDimension); ++eIter)
    {
      if (subDomain(*eIter))
      {
        const std::size_t cid = m.getContainingCell(*eIter);
        const MeshEntity localEntity = m.getLocalEntity(cid, *eIter); 
        localMatrix.clear();
        
        const typename std::map<MeshEntity, cell_integral_t>::const_iterator matIter = localMatrices.find(localEntity);
  
        if (matIter != localMatrices.end())
        {
          // Build concrete local assembly matrix
          matIter->second.evaluate(localMatrix, cid, values);
          addValues(cid, localMatrix);
        }
        else
        {
          CFD_EXCEPTION("Missing optimised local assembly matrix for requested mesh entity");
        }
      }
    }

    assemble();
    std::cout << "done." << std::endl;
  }


  void addToDiagonal(DiscreteField<dimension>& v)
  {
    assert(rowMappings == v.getRowMappings());
    matrix.addToDiagonal(v.getVectorHandle());
  }

  void zeroRow(const dof_t& dof, const double diagonal)
  {
    const int rowIndex = rowMappings.getGlobalIndex(dof);
    matrix.zeroRow(rowIndex, diagonal);
  }

  void zeroRows(const DofMap<dimension>& dofs, const double diagonal)
  { 
    std::vector<int> rowIndices = dofs.getIndices(rowMappings);
    matrix.zeroRows(rowIndices.size(), &rowIndices[0], diagonal);
  }
  
  void zero()
  {
    matrix.zero();
  }

  DiscreteField<dimension> getLumpedDiagonal() const
  {
    assert(rowMappings == colMappings);
    return DiscreteField<dimension>(rowMappings, matrix.getLumpedDiagonal());
  }

  void scaleDiagonal(const DiscreteField<dimension>& s)
  {
    assert(rowMappings == colMappings);
    assert(rowMappings == s.getRowMappings());
    matrix.scaleDiagonal(s.getVectorHandle());
  }

  void assemble()
  {
    matrix.assemble();
  }

  DiscreteField<dimension> operator*(const DiscreteField<dimension>& v) const
  {
    assert(colMappings == v.getRowMappings());
    return DiscreteField<dimension>(rowMappings, matrix*v.getVectorHandle());
  }

  DiscreteOperator<dimension> operator*(const DiscreteOperator<dimension>& b) const
  {
    assert(colMappings == b.getRowMappings());
    return DiscreteOperator<dimension>(rowMappings, b.getColMappings(), matrix*b.getMatrixHandle());
  }

  DiscreteField<dimension> trans_mult(const DiscreteField<dimension>& v) const
  {
    assert(rowMappings == v.getRowMappings());
    return DiscreteField<dimension>(colMappings, matrix.trans_mult(v.getVectorHandle()));
  }

  DiscreteOperator<dimension> trans_mult(const DiscreteOperator<dimension>& b) const
  {
    assert(rowMappings == b.getRowMappings());
    return DiscreteOperator<dimension>(colMappings, b.getColMappings(), matrix.trans_mult(b.getMatrixHandle()));
  }

  PETScMatrix& getMatrixHandle()
  {
    return matrix;
  }

  const PETScMatrix& getMatrixHandle() const
  {
    return matrix;
  }
};

}

#endif
