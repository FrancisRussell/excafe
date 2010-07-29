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
#include "capture/forms/form_evaluator.hpp"
#include "capture/assembly/assembly_helper.hpp"
#include "capture/assembly/scalar_placeholder.hpp"
#include "capture/assembly/scalar_placeholder_evaluator.hpp"

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
  typedef detail::ScalarPlaceholder::expression_t                        expression_t;
  typedef detail::ScalarPlaceholder::optimised_expression_t              optimised_expression_t;
  typedef detail::LocalAssemblyMatrix<dimension, expression_t>           local_matrix_t;
  typedef detail::LocalAssemblyMatrix<dimension, optimised_expression_t> opt_local_matrix_t;
  typedef detail::LocalAssemblyMatrix<dimension, double>                 evaluated_local_matrix_t;

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

  void addTermGeneral(const Scenario<dimension>& scenario, const detail::ExpressionValues<dimension>& values,
    const forms::BilinearFormIntegralSum::const_iterator sumBegin, 
    const forms::BilinearFormIntegralSum::const_iterator sumEnd,
    const MeshFunction<bool>& subDomain)
  {
    using namespace cfd::forms;

    const std::set<const finite_element_t*> trialElements(colMappings.getFiniteElements());
    const std::set<const finite_element_t*> testElements(rowMappings.getFiniteElements());

    typedef detail::LocalAssemblyMatrix<dimension, double> local_matrix_t;
    local_matrix_t localMatrix(testElements, trialElements);

    typedef std::pair<const finite_element_t*, const finite_element_t*> element_pair;
    typedef std::pair< FormEvaluator<dimension>, FormEvaluator<dimension> > evaluator_pair;

    std::map< element_pair, std::vector<evaluator_pair> > evaluators;

    for(BilinearFormIntegralSum::const_iterator formIter = sumBegin; formIter!=sumEnd; ++formIter)
    {
      // Find trial
      BasisFinder<dimension> trialFinder(scenario);
      formIter->getTrialField()->accept(trialFinder);

      const finite_element_t* const trialBasis = trialFinder.getBasis();
      assert(trialBasis != NULL);
      assert(trialElements.find(trialBasis) != trialElements.end());

      // Find test
      BasisFinder<dimension> testFinder(scenario);
      formIter->getTestField()->accept(testFinder);

      const finite_element_t* const testBasis = testFinder.getBasis();
      assert(testBasis != NULL);
      assert(testElements.find(testBasis) != testElements.end());

      const FormEvaluator<dimension> trialEvaluator(trialBasis->getCell(), scenario, values, formIter->getTrialField());
      const FormEvaluator<dimension> testEvaluator(testBasis->getCell(), scenario, values, formIter->getTestField());
      evaluators[std::make_pair(trialBasis, testBasis)].push_back(std::make_pair(trialEvaluator, testEvaluator));
    }

    const Mesh<dimension>& m = rowMappings.getMesh();
    const std::size_t degree = 5;
    boost::array<std::size_t, dimension> degrees;
    std::fill(degrees.begin(), degrees.end(), degree);
    const QuadraturePoints<dimension> quadrature = m.getReferenceCell()->getQuadrature(degrees);

    const std::size_t entityDimension = subDomain.getDimension();

    for(typename Mesh<dimension>::global_iterator eIter(m.global_begin(entityDimension)); eIter != m.global_end(entityDimension); ++eIter)
    {
      if (subDomain(*eIter))
      {
        //FIXME: Assumes constant jacobian
        const std::size_t cid = m.getContainingCell(*eIter);
        const CellVertices<dimension> vertices(m.getCoordinates(cid));
        const MeshEntity localEntity = m.getLocalEntity(cid, *eIter); 
        const double jacobian = m.getReferenceCell()->getJacobian(vertices, localEntity, vertex_type(0.0, 0.0));
        localMatrix.clear();
  
        for(typename std::map< element_pair, std::vector<evaluator_pair> >::const_iterator evaluatorIter=evaluators.begin(); evaluatorIter!=evaluators.end(); ++evaluatorIter)
        {
          const finite_element_t* const trialFunction = evaluatorIter->first.first;
          const finite_element_t* const testFunction = evaluatorIter->first.second;
  
          const unsigned trialSpaceDimension = trialFunction->spaceDimension();
          const unsigned testSpaceDimension = testFunction->spaceDimension();
  
          std::vector<int> trialIndices(trialSpaceDimension);
          std::vector<int> testIndices(testSpaceDimension);
  
          std::vector< Tensor<dimension> > trialValues(trialSpaceDimension);
          std::vector< Tensor<dimension> > testValues(testSpaceDimension);
  
          for(typename std::vector<evaluator_pair>::const_iterator bFormIter(evaluatorIter->second.begin()); bFormIter!=evaluatorIter->second.end(); ++bFormIter)
          {
  
            for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
              trialIndices[trial] = localMatrix.getTrialOffset(*trialFunction, trial);
  
            for(unsigned test=0; test<testSpaceDimension; ++test)
              testIndices[test] = localMatrix.getTestOffset(*testFunction, test);

            for(typename QuadraturePoints<dimension>::iterator quadIter(quadrature.begin(localEntity)); quadIter!=quadrature.end(localEntity); ++quadIter)
            {
              for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
                trialValues[trial] = bFormIter->first.evaluate(vertices, localEntity, quadIter->first, Dof<dimension>(trialFunction, cid, trial));
  
              for(unsigned test=0; test<testSpaceDimension; ++test)
                testValues[test] = bFormIter->second.evaluate(vertices, localEntity, quadIter->first, Dof<dimension>(testFunction, cid, test));
  
              for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
              {
                for(unsigned test=0; test<testSpaceDimension; ++test)
                {
                  localMatrix(testIndices[test], trialIndices[trial]) += trialValues[trial].colon_product(testValues[test]) *
                    quadIter->second * jacobian;
                }
              }
            }
          }
        }
        addValues(cid, localMatrix);
      }
    }
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

  void assembleForms(Scenario<dimension>& scenario, detail::ExpressionValues<dimension>& values, const forms::BilinearFormIntegralSum& expr)
  {
    const Mesh<dimension>& m(scenario.getMesh());
    const MeshFunction<bool> allCells(dimension, true);
      addTermGeneral(scenario, values, expr.begin_dx(), expr.end_dx(), allCells);

    const MeshFunction<bool> boundaryFunction = m.getBoundaryFunction();
    addTermGeneral(scenario, values, expr.begin_ds(), expr.end_ds(), boundaryFunction);

    //FIXME: perform internal boundary integrals
    assemble();
  }

  void assembleFromOptimisedLocalMatrix(const Scenario<dimension>& scenario, 
    const detail::ExpressionValues<dimension>& values,
    const std::map<MeshEntity, opt_local_matrix_t> localMatrices,
    const MeshFunction<bool>& subDomain)
  {
    using namespace detail;

    // Build set of placeholders
    PolynomialVariableCollector<expression_t::optimised_t> collector;

    typedef std::pair<MeshEntity, opt_local_matrix_t> entity_matrix_pair_t;
    BOOST_FOREACH(const entity_matrix_pair_t& mapping, localMatrices)
    {
      collector = std::for_each(mapping.second.begin(), mapping.second.end(), collector);
    }

    const std::set<ScalarPlaceholder> placeholders(collector.getVariables());
    const std::set<const finite_element_t*> trialElements(colMappings.getFiniteElements());
    const std::set<const finite_element_t*> testElements(rowMappings.getFiniteElements());
    const Mesh<dimension>& m = scenario.getMesh();
    const std::size_t entityDimension = subDomain.getDimension();
    local_matrix_t localMatrix(testElements, trialElements);

    std::cout << "Starting assembly...." << std::flush;
    for(typename Mesh<dimension>::global_iterator eIter(m.global_begin(entityDimension)); eIter != m.global_end(entityDimension); ++eIter)
    {
      if (subDomain(*eIter))
      {
        const std::size_t cid = m.getContainingCell(*eIter);
        const MeshEntity localEntity = m.getLocalEntity(cid, *eIter); 
        localMatrix.clear();
        
        const typename std::map<MeshEntity, opt_local_matrix_t>::const_iterator matIter = localMatrices.find(localEntity);
  
        if (matIter != localMatrices.end())
        {
          // Find placeholder values
          typedef optimised_expression_t::value_map value_map;
          const value_map placeholderValues(evaluatePlaceholders<value_map>(scenario, values, cid, placeholders));

          // Build concrete local assembly matrix
          const PolynomialEvaluator<optimised_expression_t> evaluator(placeholderValues);
          const evaluated_local_matrix_t concreteLocalMatrix(matIter->second.transform(evaluator));
          addValues(cid, concreteLocalMatrix);
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

  void extractSubmatrix(DiscreteOperator& s) const
  {
    const std::vector<int> rowIndices = s.rowMappings.getIndices(rowMappings);
    const std::vector<int> colIndices = s.colMappings.getIndices(colMappings);
    matrix.extractSubmatrix(s.matrix, rowIndices.size(), colIndices.size(), &rowIndices[0], &colIndices[0]);
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
