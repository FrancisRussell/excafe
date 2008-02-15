#ifndef SIMPLE_CFD_NUMERIC_MATRIX_HPP
#define SIMPLE_CFD_NUMERIC_MATRIX_HPP

#include <cassert>
#include "petsc.h"
#include "petscmat.h"
#include "petscis.h"

namespace cfd
{

class PETScMatrix
{
private:
  Mat m;

  void checkError(const PetscErrorCode ierr) const
  {
    assert(ierr == 0);
  }

public:
  PETScMatrix(const unsigned rows, const unsigned cols)
  {
    PetscErrorCode ierr;
    ierr = MatCreateSeqAIJ(PETSC_COMM_SELF, rows, cols, 0, PETSC_NULL, &m);
    checkError(ierr);
  }

  void addValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block)
  {
    const PetscErrorCode ierr = MatSetValues(m, rows, rowIndices, cols, colIndices, block, ADD_VALUES);
    checkError(ierr);
  }

  void setValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, const double* block)
  {
    const PetscErrorCode ierr = MatSetValues(m, rows, rowIndices, cols, colIndices, block, INSERT_VALUES);
    checkError(ierr);
  }

  void getValues(const unsigned rows, const unsigned cols, const int* rowIndices, const int* colIndices, double* block) const
  {
    const PetscErrorCode ierr = MatGetValues(m, rows, rowIndices, cols, colIndices, block);
    checkError(ierr);
  }

  void zeroRow(const int row)
  {
    zeroRows(&row, 1);
  }

  void zeroRows(const int* rows, const unsigned rowCount)
  {
    IS indexSet;
    PetscErrorCode ierr = ISCreateGeneral(PETSC_COMM_SELF, rowCount, rows, &indexSet);
    checkError(ierr);
    
    ierr = MatZeroRowsIS(m, indexSet, 0.0);
    checkError(ierr);

    ierr = ISDestroy(indexSet);
    checkError(ierr);
  }

  void assemble()
  {
    MatAssemblyBegin(m, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m, MAT_FINAL_ASSEMBLY);
  }

  ~PETScMatrix()
  {
    const PetscErrorCode ierr = MatDestroy(m);
    checkError(ierr);
  }
};

}

#endif
