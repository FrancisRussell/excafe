#ifndef EXCAFE_NUMERIC_MATRIX_HPP
#define EXCAFE_NUMERIC_MATRIX_HPP

#include "numeric/sparsity_pattern.hpp"
#include "numeric/vector.hpp"
#include "petsc.h"
#include "petscmat.h"
#include "petscis.h"

namespace excafe
{

class PETScMatrix
{
private:
  Mat m;

  PETScMatrix(const Mat& _m);
  void checkError(const PetscErrorCode ierr) const;

public:
  PETScMatrix(const PETScMatrix& p);
  PETScMatrix(const unsigned rows, const unsigned cols);
  PETScMatrix(const unsigned rows, const unsigned cols, const unsigned nz);
  PETScMatrix(const SparsityPattern& pattern);
  PETScMatrix& operator=(const PETScMatrix& r);
  PETScVector operator*(const PETScVector& v) const;
  PETScMatrix operator*(const PETScMatrix& m) const;
  PETScVector trans_mult(const PETScVector& m) const;
  PETScMatrix trans_mult(const PETScMatrix& m) const;
  std::size_t numRows() const;
  std::size_t numCols() const;
  void addValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block);
  void setValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block);
  void getValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, double* block) const;
  void zeroRow(const int row, const double diagonal);
  void zeroRows(const unsigned rowCount, const int* rows, const double diagonal);
  void addToDiagonal(const PETScVector& v);
  void scaleDiagonal(const PETScVector& s);
  void zero();
  void assemble();
  void view() const;
  PETScVector getLumpedDiagonal() const;
  Mat getPETScHandle() const;
  ~PETScMatrix();
};

}

#endif
