#ifndef SIMPLE_CFD_FE_MATRIX_HPP
#define SIMPLE_CFD_FE_MATRIX_HPP

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include "vertex.hpp"
#include "dof_map.hpp"
#include "fe_binary_function.hpp"
#include "fe_vector.hpp"
#include "mesh.hpp"
#include "numeric/matrix.hpp"
#include "numeric/sparsity_pattern.hpp"
#include "quadrature_points.hpp"
#include "forms/bilinear_form_sum.hpp"
#include "forms/basis_finder.hpp"
#include "forms/field.hpp"
#include "forms/form_evaluator.hpp"

namespace cfd
{

template<std::size_t D>
class FEMatrix
{
private:
  static const std::size_t dimension = D;
  typedef vertex<dimension> vertex_type;
  typedef FiniteElement<dimension> finite_element_t;
  typedef typename DofMap<dimension>::dof_t dof_t;

  const DofMap<dimension> rowMappings;
  const DofMap<dimension> colMappings;
  PETScMatrix matrix;

  static SparsityPattern createSparsityPattern(const DofMap<dimension>& rowMappings,
                                               const DofMap<dimension>& colMappings)
  {
    const unsigned rowDofs = rowMappings.getDegreesOfFreedomCount();
    const unsigned colDofs = colMappings.getDegreesOfFreedomCount();

    SparsityPattern pattern(rowDofs, colDofs);
    const Mesh<dimension> m(rowMappings.getMesh());

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
                colMappings.getGlobalIndexWithMissingAsNegative(dof_t(*rowElemIter,
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

  template<typename cell_type>
  void addTermGeneral(const Mesh<dimension>& m, const FEBinaryFunction<cell_type>& f, const bool boundaryIntegral)
  {
    const std::set<const finite_element_t*> trialElements(colMappings.getFiniteElements());
    const std::set<const finite_element_t*> testElements(rowMappings.getFiniteElements());
    
    const finite_element_t* const trialFunction = f.getTrialFunction();
    const finite_element_t* const testFunction = f.getTestFunction();

    assert(trialElements.find(trialFunction) != trialElements.end());
    assert(testElements.find(testFunction) != testElements.end());

    const std::size_t entityDimension = boundaryIntegral ? m.getDimension()-1 : m.getDimension();
    const unsigned testSpaceDimension = testFunction->spaceDimension();
    const unsigned trialSpaceDimension = trialFunction->spaceDimension();
    const MeshFunction<bool> boundaryFunction = m.getBoundaryFunction();

    std::vector<int> testIndices(testSpaceDimension);
    std::vector<int> trialIndices(trialSpaceDimension);
    std::vector<double> valueBlock(testSpaceDimension*trialSpaceDimension);

    const std::size_t degree = 5;
    const QuadraturePoints<dimension> quadrature = m.getReferenceCell().getQuadrature(degree);

    for(typename Mesh<dimension>::global_iterator eIter(m.global_begin(entityDimension)); eIter != m.global_end(entityDimension); ++eIter)
    {
      if (!boundaryIntegral || boundaryFunction(*eIter))
      {
        //FIXME: Assumes constant jacobian
        const std::size_t cid = m.getContainingCell(*eIter);
        const CellVertices<dimension> vertices(m.getCoordinates(cid));
        const MeshEntity localEntity = m.getLocalEntity(cid, *eIter); 
        const double jacobian = m.getReferenceCell().getJacobian(vertices, localEntity, vertex_type(0.0, 0.0));

        std::fill(valueBlock.begin(), valueBlock.end(), 0.0);

        for(unsigned test=0; test<testSpaceDimension; ++test)
          testIndices[test] = rowMappings.getGlobalIndexWithMissingAsNegative(dof_t(testFunction, cid, test));

        for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
          trialIndices[trial] = colMappings.getGlobalIndexWithMissingAsNegative(dof_t(trialFunction, cid, trial));

        for(typename QuadraturePoints<dimension>::iterator quadIter(quadrature.begin(localEntity)); quadIter!=quadrature.end(localEntity); ++quadIter)
          for(unsigned test=0; test<testSpaceDimension; ++test)
            for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
              valueBlock[test * trialSpaceDimension + trial] += f.evaluate(vertices, *eIter, localEntity, test, trial,
              quadIter->first) * quadIter->second * jacobian;

        matrix.addValues(testSpaceDimension, trialSpaceDimension, &testIndices[0], &trialIndices[0], &valueBlock[0]);
      }
    }
  }

public:
  FEMatrix(const DofMap<dimension>& _rowMappings, const DofMap<dimension>& _colMappings) :
          rowMappings(_rowMappings), colMappings(_colMappings), matrix(createSparsityPattern(rowMappings, colMappings))
  {
    // Make sure the degree-of-freedom maps are defined on the same mesh. Calculating a
    // sparsity pattern doesn't make sense otherwise.
    assert(&rowMappings.getMesh() == &colMappings.getMesh());
  }

  FEMatrix(const DofMap<dimension>& _rowMappings, const DofMap<dimension>& _colMappings, const PETScMatrix& m) :
          rowMappings(_rowMappings), colMappings(_colMappings), matrix(m)
  {
    assert(&rowMappings.getMesh() == &colMappings.getMesh());
  }

  FEMatrix(const FEMatrix& m) : rowMappings(m.rowMappings), colMappings(m.colMappings), matrix(m.matrix)
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

  FEMatrix& operator=(const FEMatrix& f)
  {
    assert(rowMappings == f.rowMappings);
    assert(colMappings == f.colMappings);
    matrix = f.matrix;
    return *this;
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

  template<typename cell_type>
  void addTerm(const Mesh<dimension>& m, const FEBinaryFunction<cell_type>& f)
  {
    addTermGeneral(m, f, false);
  }

  template<typename cell_type>
  void addBoundaryTerm(const Mesh<dimension>& m, const FEBinaryFunction<cell_type>& f)
  {
    addTermGeneral(m, f, true);
  }

  FEMatrix& operator+=(const forms::BilinearFormSum& expr)
  {
    using namespace cfd::forms;

    const std::set<const finite_element_t*> trialElements(colMappings.getFiniteElements());
    const std::set<const finite_element_t*> testElements(rowMappings.getFiniteElements());

    typedef std::pair<const finite_element_t*, const finite_element_t*> element_pair;
    typedef std::pair< FormEvaluator<dimension>, FormEvaluator<dimension> > evaluator_pair;

    std::map< element_pair, std::vector<evaluator_pair> > evaluators;

    for(BilinearFormSum::const_iterator formIter = expr.begin(); formIter!=expr.end(); ++formIter)
    {
      // Find trial
      BasisFinder<dimension> trialFinder;
      formIter->getTrialField()->accept(trialFinder);

      const finite_element_t* const trialBasis = trialFinder.getBasis();
      assert(trialBasis != NULL);
      assert(trialElements.find(trialBasis) != trialElements.end());

      // Find test
      BasisFinder<dimension> testFinder;
      formIter->getTestField()->accept(testFinder);

      const finite_element_t* const testBasis = testFinder.getBasis();
      assert(testBasis != NULL);
      assert(testElements.find(testBasis) != testElements.end());

      const FormEvaluator<dimension> trialEvaluator(formIter->getTrialField(), trialBasis->getCell());
      const FormEvaluator<dimension> testEvaluator(formIter->getTestField(), testBasis->getCell());
      evaluators[std::make_pair(trialBasis, testBasis)].push_back(std::make_pair(trialEvaluator, testEvaluator));
    }

    return *this;
  }

  void addToDiagonal(FEVector<dimension>& v)
  {
    assert(rowMappings == v.getRowMappings());
    matrix.addToDiagonal(v.getVectorHandle());
  }

  void zeroRow(const dof_t& dof, const double diagonal)
  {
    const int rowIndex = rowMappings.getGlobalIndex(dof);
    matrix.zeroRow(rowIndex, diagonal);
  }
  
  void zero()
  {
    matrix.zero();
  }

  void extractSubmatrix(FEMatrix& s) const
  {
    const std::vector<int> rowIndices = s.rowMappings.getIndices(rowMappings);
    const std::vector<int> colIndices = s.colMappings.getIndices(colMappings);
    matrix.extractSubmatrix(s.matrix, rowIndices.size(), colIndices.size(), &rowIndices[0], &colIndices[0]);
  }

  FEVector<dimension> getLumpedDiagonal() const
  {
    assert(rowMappings == colMappings);
    return FEVector<dimension>(rowMappings, matrix.getLumpedDiagonal());
  }

  void scaleDiagonal(const FEVector<dimension>& s)
  {
    assert(rowMappings == colMappings);
    assert(rowMappings == s.getRowMappings());
    matrix.scaleDiagonal(s.getVectorHandle());
  }

  void assemble()
  {
    matrix.assemble();
  }

  FEVector<dimension> operator*(const FEVector<dimension>& v) const
  {
    assert(colMappings == v.getRowMappings());
    return FEVector<dimension>(rowMappings, matrix*v.getVectorHandle());
  }

  FEMatrix<dimension> operator*(const FEMatrix<dimension>& b) const
  {
    assert(colMappings == b.getRowMappings());
    return FEMatrix<dimension>(rowMappings, b.getColMappings(), matrix*b.getMatrixHandle());
  }

  FEVector<dimension> trans_mult(const FEVector<dimension>& v) const
  {
    assert(rowMappings == v.getRowMappings());
    return FEVector<dimension>(colMappings, matrix.trans_mult(v.getVectorHandle()));
  }

  FEMatrix<dimension> trans_mult(const FEMatrix<dimension>& b) const
  {
    assert(rowMappings == b.getRowMappings());
    return FEMatrix<dimension>(colMappings, b.getColMappings(), matrix.trans_mult(b.getMatrixHandle()));
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
