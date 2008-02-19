#ifndef SIMPLE_CFD_NUMERIC_VECTOR_HPP
#define SIMPLE_CFD_NUMERIC_VECTOR_HPP

#include "petsc.h"
#include "petscvec.h"

namespace cfd
{

class PETScVector
{
private:
  Vec v;

  void checkError(const PetscErrorCode ierr) const
  {
    assert(ierr == 0);
  }

public:
  PETScVector(const unsigned rows)
  {
    const PetscErrorCode ierr = VecCreateSeq(PETSC_COMM_SELF, rows, &v);
    checkError(ierr);
  }

  void addValues(const unsigned numValues, const int* indices, const double* values)
  {
    const PetscErrorCode ierr = VecSetValues(v, numValues, indices, values, ADD_VALUES);
    checkError(ierr);
  }

  void setValues(const unsigned numValues, const int* indices, const double* values)
  {
    const PetscErrorCode ierr = VecSetValues(v, numValues, indices, values, INSERT_VALUES);
    checkError(ierr);
  }

  void getValues(const unsigned numValues, const int* indices, double* values) const
  {
    const PetscErrorCode ierr = VecGetValues(v, numValues, indices, values);
    checkError(ierr);
  }

  void assemble()
  {
    VecAssemblyBegin(v);
    VecAssemblyEnd(v);
  }

  Vec getPETScHandle()
  {
    return v;
  }

  ~PETScVector()
  {
    const PetscErrorCode ierr = VecDestroy(v);
    checkError(ierr);
  }
};

}

#endif
