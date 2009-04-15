#include "simple_cfd/numeric/matrix.hpp"
#include "simple_cfd/numeric/vector.hpp"
#include "simple_cfd/numeric/sparsity_pattern.hpp"
#include <vector>
#include <cassert>
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
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

PETScVector PETScMatrix::operator*(const PETScVector& v) const
{
  PETScVector result(numRows());
  const PetscErrorCode ierr = MatMult(m, v.getPETScHandle(), result.getPETScHandle());
  checkError(ierr);
  return result;
}

PETScMatrix PETScMatrix::operator*(const PETScMatrix& b) const
{
  Mat c;
  const PetscErrorCode ierr = MatMatMult(m, b.getPETScHandle(), MAT_INITIAL_MATRIX, 1.0, &c);
  checkError(ierr);
  return PETScMatrix(c);
}

PETScVector PETScMatrix::trans_mult(const PETScVector& v) const
{
  PETScVector result(numCols());
  const PetscErrorCode ierr = MatMultTranspose(m, v.getPETScHandle(), result.getPETScHandle());
  checkError(ierr);
  return result;
}

PETScMatrix PETScMatrix::trans_mult(const PETScMatrix& b) const
{
  Mat c;
  const PetscErrorCode ierr = MatMatMultTranspose(m, b.getPETScHandle(), MAT_INITIAL_MATRIX, 1.0, &c);
  checkError(ierr);
  return PETScMatrix(c);
}

std::size_t PETScMatrix::numRows() const
{
  PetscInt rows;
  const PetscErrorCode ierr = MatGetSize(m, &rows, PETSC_NULL);
  checkError(ierr);

  return rows;
}

std::size_t PETScMatrix::numCols() const
{
  PetscInt cols;
  const PetscErrorCode ierr = MatGetSize(m, PETSC_NULL, &cols);
  checkError(ierr);

  return cols;
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

void PETScMatrix::addToDiagonal(const PETScVector& v)
{
  assert(numRows() == numCols());
  assert(numRows() == v.numRows());

  const PetscErrorCode ierr = MatDiagonalSet(m, v.getPETScHandle(), ADD_VALUES);
  checkError(ierr);
}

void PETScMatrix::zero()
{
  const PetscErrorCode ierr = MatZeroEntries(m);
  checkError(ierr);
}

void PETScMatrix::assemble()
{
  MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
}

Mat PETScMatrix::getPETScHandle() const
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

PETScVector PETScMatrix::getLumpedDiagonal() const
{
  assert(numRows() == numCols());

  PETScVector allOnes(numCols());
  allOnes = 1.0;
  allOnes.assemble();

  return *this * allOnes;
}

void PETScMatrix::extractSubmatrix(PETScMatrix& dest, const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices) const
{
  const PetscErrorCode ierr = MatGetSubMatrixRaw(m, rows, rowIndices, cols, colIndices, PETSC_DECIDE, MAT_INITIAL_MATRIX, &dest.m);
  checkError(ierr);
}

PETScMatrix::~PETScMatrix()
{
  const PetscErrorCode ierr = MatDestroy(m);
  checkError(ierr);
}

}

