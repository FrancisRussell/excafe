#include "simple_cfd/numeric/matrix.hpp"
#include "simple_cfd/numeric/sparsity_pattern.hpp"
#include <vector>
#include <cassert>
#include "petsc.h"
#include "petscmat.h"
#include "petscis.h"

namespace cfd
{

PETScMatrix::PETScMatrix(const Mat& _m) : m(_m)
{
}

PETScMatrix::PETScMatrix(const PETScMatrix& p)
{
  const PetscErrorCode ierr = MatConvert(p.m, MATSAME, MAT_INITIAL_MATRIX, &m);
  checkError(ierr);
}

void PETScMatrix::checkError(const PetscErrorCode ierr) const
{
  assert(ierr == 0);
}

PETScMatrix::PETScMatrix(const unsigned rows, const unsigned cols)
{
  const PetscErrorCode ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, PETSC_DEFAULT, PETSC_NULL, &m);
  checkError(ierr);
}

PETScMatrix::PETScMatrix(const unsigned rows, const unsigned cols, const unsigned nz)
{
  const PetscErrorCode ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, nz, PETSC_NULL, &m);
  checkError(ierr);
}

PETScMatrix::PETScMatrix(const SparsityPattern& pattern)
{
  const unsigned rows = pattern.getRows();
  const unsigned cols = pattern.getCols();
  const std::vector<int> nnz(pattern.getNonZerosPerRow());
  const PetscErrorCode ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, PETSC_DEFAULT, &nnz[0], &m);
  checkError(ierr);
}

PETScMatrix& PETScMatrix::operator=(const PETScMatrix& r)
{
  const PetscErrorCode ierr = MatCopy(r.m, m, DIFFERENT_NONZERO_PATTERN);
  checkError(ierr);
  return *this;
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

void PETScMatrix::view() const
{
  PetscErrorCode ierr = PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
  checkError(ierr);

  ierr = MatView(m, PETSC_VIEWER_STDOUT_WORLD);
  checkError(ierr);
}

PETScMatrix PETScMatrix::extractSubmatrix(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices) const
{
  Mat submatrix;
  const PetscErrorCode ierr = MatGetSubMatrixRaw(m, rows, rowIndices, cols, colIndices, PETSC_DECIDE, MAT_INITIAL_MATRIX, &submatrix);
  checkError(ierr);
  return PETScMatrix(submatrix);
}

PETScMatrix::~PETScMatrix()
{
  const PetscErrorCode ierr = MatDestroy(m);
  checkError(ierr);
}

}

