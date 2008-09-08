#ifndef SIMPLE_CFD_NUMERIC_MATRIX_HPP
#define SIMPLE_CFD_NUMERIC_MATRIX_HPP

#include "numeric/sparsity_pattern.hpp"
#include "petsc.h"
#include "petscmat.h"
#include "petscis.h"

namespace cfd
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
  void addValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block);
  void setValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block);
  void getValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, double* block) const;
  void zeroRow(const int row, const double diagonal);
  void zeroRows(const int* rows, const unsigned rowCount, const double diagonal);
  void assemble();
  PETScMatrix extractSubmatrix(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices) const;
  void view() const;
  Mat getPETScHandle();
  ~PETScMatrix();
};

}

#endif
