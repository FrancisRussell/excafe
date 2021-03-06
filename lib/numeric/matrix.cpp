#include <boost/static_assert.hpp>
#include "excafe/numeric/matrix.hpp"
#include "excafe/numeric/vector.hpp"
#include "excafe/numeric/sparsity_pattern.hpp"
#include <vector>
#include <cassert>
#include "petsc.h"
#include "petscmat.h"
#include "petscvec.h"
#include "petscis.h"

namespace excafe
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
  const PetscErrorCode ierr = MatMatTransposeMult(m, b.getPETScHandle(), MAT_INITIAL_MATRIX, 1.0, &c);
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
  zeroRows(1, &row, diagonal);
}

// Handle PETSc API change between 3.1 and 3.2.

template<typename ZeroFunction>
PetscErrorCode ZeroRowsWrapper(ZeroFunction func, Mat& mat, const unsigned numRows, const int* const rows, const double diag)
{
  BOOST_STATIC_ASSERT(sizeof(ZeroFunction) == 0);
  return 0;
}

template<>
PetscErrorCode ZeroRowsWrapper(PetscErrorCode (*func)(Mat, PetscInt, const PetscInt*, PetscScalar, Vec, Vec),
    Mat& mat, const unsigned numRows, const int* const rows, const double diag)
{
  return func(mat, numRows, rows, diag, PETSC_NULL, PETSC_NULL);
}

template<>
PetscErrorCode ZeroRowsWrapper(PetscErrorCode (*func)(Mat, PetscInt, const PetscInt*, PetscScalar),
    Mat& mat, const unsigned numRows, const int* const rows, const double diag)
{
  return func(mat, numRows, rows, diag);
}

void PETScMatrix::zeroRows(const unsigned rowCount, const int* rows, const double diagonal)
{
  const PetscErrorCode ierr = ZeroRowsWrapper(&MatZeroRows, m, rowCount, rows, diagonal);
  checkError(ierr);
}

void PETScMatrix::addToDiagonal(const PETScVector& v)
{
  assert(numRows() == numCols());
  assert(numRows() == v.numRows());

  const PetscErrorCode ierr = MatDiagonalSet(m, v.getPETScHandle(), ADD_VALUES);
  checkError(ierr);
}

void PETScMatrix::scaleDiagonal(const PETScVector& s)
{
  assert(numRows() == numCols());
  assert(numRows() == s.numRows());

  PETScVector diagonal(s.numRows());

  PetscErrorCode ierr = MatGetDiagonal(m, diagonal.getPETScHandle());
  checkError(ierr);

  ierr = VecPointwiseMult(diagonal.getPETScHandle(), diagonal.getPETScHandle(), s.getPETScHandle());
  checkError(ierr);

  ierr = MatDiagonalSet(m, diagonal.getPETScHandle(), INSERT_VALUES);
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


// Handle PETSc API change between 3.1 and 3.2.

template<typename DestroyFunction>
PetscErrorCode MatDestroyWrapper(DestroyFunction func, Mat& m)
{
  BOOST_STATIC_ASSERT(sizeof(DestroyFunction) == 0);
  return 0;
}

template<>
PetscErrorCode MatDestroyWrapper(PetscErrorCode (*func)(Mat*), Mat& m)
{
  return func(&m);
}

template<>
PetscErrorCode MatDestroyWrapper(PetscErrorCode (*func)(Mat), Mat& m)
{
  return func(m);
}

PETScMatrix::~PETScMatrix()
{
  const PetscErrorCode ierr = MatDestroyWrapper(&MatDestroy, m);
  checkError(ierr);
}

}

