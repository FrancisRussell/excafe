#include "simple_cfd/numeric/matrix.hpp"
#include "simple_cfd/numeric/sparsity_pattern.hpp"
#include <cassert>
#include "petsc.h"
#include "petscmat.h"
#include "petscis.h"

namespace cfd
{

void PETScMatrix::checkError(const PetscErrorCode ierr) const
{
  assert(ierr == 0);
}

PETScMatrix::PETScMatrix(const unsigned rows, const unsigned cols)
{
  PetscErrorCode ierr;
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, PETSC_DEFAULT, PETSC_NULL, &m);
  checkError(ierr);
}

PETScMatrix::PETScMatrix(const unsigned rows, const unsigned cols, const unsigned nz)
{
  PetscErrorCode ierr;
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, nz, PETSC_NULL, &m);
  checkError(ierr);
}

PETScMatrix::PETScMatrix(const SparsityPattern& pattern)
{
  const unsigned rows = pattern.getRows();
  const unsigned cols = pattern.getCols();
  const std::vector<int> nnz(pattern.getNonZerosPerRow());
  PetscErrorCode ierr;
  ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, PETSC_DEFAULT, &nnz[0], &m);
  checkError(ierr);
}

void PETScMatrix::addValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block)
{
  const PetscErrorCode ierr = MatSetValues(m, rows, rowIndices, cols, colIndices, block, ADD_VALUES);
  checkError(ierr);
}

void PETScMatrix::setValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block)
{
  const PetscErrorCode ierr = MatSetValues(m, rows, rowIndices, cols, colIndices, block, INSERT_VALUES);
  checkError(ierr);
}

void PETScMatrix::getValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, double* block) const
{
  const PetscErrorCode ierr = MatGetValues(m, rows, rowIndices, cols, colIndices, block);
  checkError(ierr);
}

void PETScMatrix::zeroRow(const int row, const double diagonal)
{
  zeroRows(&row, 1, diagonal);
}

void PETScMatrix::zeroRows(const int* rows, const unsigned rowCount, const double diagonal)
{
  IS indexSet;
  PetscErrorCode ierr = ISCreateGeneral(PETSC_COMM_SELF, rowCount, rows, &indexSet);
  checkError(ierr);
    
  ierr = MatZeroRowsIS(m, indexSet, diagonal);
  checkError(ierr);

  ierr = ISDestroy(indexSet);
  checkError(ierr);
}

void PETScMatrix::assemble()
{
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
}

Mat PETScMatrix::getPETScHandle()
{
  return m;
}

PETScMatrix::~PETScMatrix()
{
  const PetscErrorCode ierr = MatDestroy(m);
  checkError(ierr);
}

}

