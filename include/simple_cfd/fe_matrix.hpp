#ifndef SIMPLE_CFD_FE_MATRIX_HPP
#define SIMPLE_CFD_FE_MATRIX_HPP

#include "numeric/matrix.hpp"
#include "numeric/sparsity_pattern.hpp"
#include "dof_map.hpp"
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

  void zeroRow(const dof_t& dof, const double diagonal)
  {
    const int rowIndex = rowMappings.getGlobalIndex(dof);
    matrix.zeroRow(rowIndex, diagonal);
  }

  FEMatrix extractSubmatrix(const finite_element_t* rowElement, const finite_element_t* colElement) const
  {
    const dof_map<cell_type> newRowMappings = rowMappings.extractDofs(rowElement);
    const dof_map<cell_type> newColMappings = colMappings.extractDofs(colElement);

    std::vector<int> rowIndices = newRowMappings.getIndices(rowMappings);
    std::vector<int> colIndices = newColMappings.getIndices(colMappings);

    PETScMatrix subMatrix = matrix.extractSubmatrix(rowIndices.size(), colIndices.size(), &rowIndices[0], &colIndices[0]);

    return FEMatrix(newRowMappings, newColMappings, subMatrix);
  }

  void assemble()
  {
    matrix.assemble();
  }

  PETScMatrix& getMatrixHandle()
  {
    return matrix;
  }
};

}

#endif
