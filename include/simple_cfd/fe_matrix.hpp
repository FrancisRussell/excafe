#ifndef SIMPLE_CFD_FE_MATRIX_HPP
#define SIMPLE_CFD_FE_MATRIX_HPP

#include "fe_vector.hpp"
#include "numeric/matrix.hpp"
#include "numeric/sparsity_pattern.hpp"
#include "dof_map.hpp"
#include "fe_binary_function.hpp"
#include "mesh.hpp"
#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <boost/tuple/tuple.hpp>

namespace cfd
{

template<typename C>
class FEMatrix
{
private:
  typedef C cell_type;
  typedef typename  cell_type::vertex_type vertex_type;
  typedef finite_element<cell_type> finite_element_t;
  typedef typename dof_map<cell_type>::dof_t dof_t;

  const dof_map<cell_type> rowMappings;
  const dof_map<cell_type> colMappings;
  PETScMatrix matrix;

  static SparsityPattern createSparsityPattern(const dof_map<cell_type>& rowMappings,
                                               const dof_map<cell_type>& colMappings)
  {
    const unsigned rowDofs = rowMappings.getDegreesOfFreedomCount();
    const unsigned colDofs = colMappings.getDegreesOfFreedomCount();

    SparsityPattern pattern(rowDofs, colDofs);
    const std::map<cell_id, cell_type> cells(rowMappings.getMesh().getCells());

    const std::set<const finite_element_t*> rowElements(rowMappings.getFiniteElements());
    const std::set<const finite_element_t*> colElements(colMappings.getFiniteElements());

    for(typename std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter!=cells.end(); ++cellIter)
      for(typename std::set<const finite_element_t*>::const_iterator rowElemIter(rowElements.begin()); rowElemIter!=rowElements.end(); ++rowElemIter)
        for(typename std::set<const finite_element_t*>::const_iterator colElemIter(colElements.begin()); colElemIter!=colElements.end(); ++colElemIter)
          for(unsigned rowDof=0; rowDof < (*rowElemIter)->space_dimension(); ++rowDof)
            for(unsigned colDof=0; colDof < (*colElemIter)->space_dimension(); ++colDof)
              pattern.insert(rowMappings.getGlobalIndex(boost::make_tuple(*rowElemIter, cellIter->first, rowDof)), 
                             colMappings.getGlobalIndex(boost::make_tuple(*colElemIter, cellIter->first, colDof)));

    return pattern;
  }

public:
  FEMatrix(const dof_map<cell_type>& _rowMappings, const dof_map<cell_type>& _colMappings) :
          rowMappings(_rowMappings), colMappings(_colMappings), matrix(createSparsityPattern(rowMappings, colMappings))
  {
    // Make sure the degree-of-freedom maps are defined on the same mesh. Calculating a
    // sparsity pattern doesn't make sense otherwise.
    assert(&rowMappings.getMesh() == &colMappings.getMesh());
  }

  FEMatrix(const dof_map<cell_type>& _rowMappings, const dof_map<cell_type>& _colMappings, const PETScMatrix& m) :
          rowMappings(_rowMappings), colMappings(_colMappings), matrix(m)
  {
    assert(&rowMappings.getMesh() == &colMappings.getMesh());
  }

  FEMatrix(const FEMatrix& m) : rowMappings(m.rowMappings), colMappings(m.colMappings), matrix(m.matrix)
  {
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

  void addTerm(const mesh<cell_type>& m, const FEBinaryFunction<cell_type>& f)
  {
    const std::set<const finite_element_t*> trialElements(colMappings.getFiniteElements());
    const std::set<const finite_element_t*> testElements(rowMappings.getFiniteElements());
    
    const finite_element_t* const trialFunction = f.getTrialFunction();
    const finite_element_t* const testFunction = f.getTestFunction();

    const std::map<cell_id, cell_type> cells(m.getCells());

    const unsigned testSpaceDimension = testFunction->space_dimension();
    const unsigned trialSpaceDimension = trialFunction->space_dimension();

    std::vector<int> testIndices(testSpaceDimension);
    std::vector<int> trialIndices(trialSpaceDimension);
    std::vector<double> valueBlock(testSpaceDimension*trialSpaceDimension);

    for(typename std::map<cell_id, cell_type>::const_iterator cellIter(cells.begin()); cellIter != cells.end(); ++cellIter)
    {
      const std::map<vertex_type, double> quadrature(cellIter->second.getQuadrature(m.getGeometry()));
      std::fill(valueBlock.begin(), valueBlock.end(), 0.0);

      for(unsigned test=0; test<testSpaceDimension; ++test)
        testIndices[test] = rowMappings.getGlobalIndex(boost::make_tuple(testFunction, cellIter->first, test));

      for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
        trialIndices[trial] = colMappings.getGlobalIndex(boost::make_tuple(trialFunction, cellIter->first, trial));

      for(typename std::map<vertex_type, double>::const_iterator quadIter(quadrature.begin()); quadIter != quadrature.end(); ++quadIter)
        for(unsigned test=0; test<testSpaceDimension; ++test)
          for(unsigned trial=0; trial<trialSpaceDimension; ++trial)
            valueBlock[test * trialSpaceDimension + trial] += quadIter->second * f.evaluate(*cellIter, test, trial, quadIter->first);

      matrix.addValues(testSpaceDimension, trialSpaceDimension, &testIndices[0], &trialIndices[0], &valueBlock[0]);
    }
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

  void assemble()
  {
    matrix.assemble();
  }

  FEVector<cell_type> operator*(FEVector<cell_type>& v) const
  {
    return FEVector<cell_type>(rowMappings, matrix*v.getVectorHandle());
  }

  PETScMatrix& getMatrixHandle()
  {
    return matrix;
  }
};

}

#endif
