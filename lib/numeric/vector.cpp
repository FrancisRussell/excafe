#include "simple_cfd/numeric/vector.hpp"
#include <cassert>
#include "petsc.h"
#include "petscvec.h"

namespace cfd
{

void PETScVector::checkError(const PetscErrorCode ierr) const
{
  assert(ierr == 0);
}

PETScVector::PETScVector(const unsigned rows)
{
  const PetscErrorCode ierr = VecCreateSeq(PETSC_COMM_SELF, rows, &v);
  checkError(ierr);
}

void PETScVector::addValues(const unsigned numValues, const int* indices, const double* values)
{
  const PetscErrorCode ierr = VecSetValues(v, numValues, indices, values, ADD_VALUES);
  checkError(ierr);
}

void PETScVector::setValues(const unsigned numValues, const int* indices, const double* values)
{
  const PetscErrorCode ierr = VecSetValues(v, numValues, indices, values, INSERT_VALUES);
  checkError(ierr);
}

void PETScVector::getValues(const unsigned numValues, const int* indices, double* values) const
{
  const PetscErrorCode ierr = VecGetValues(v, numValues, indices, values);
  checkError(ierr);
}

void PETScVector::assemble()
{
  VecAssemblyBegin(v);
  VecAssemblyEnd(v);
}

Vec PETScVector::getPETScHandle()
{
  return v;
}

PETScVector::~PETScVector()
{
  const PetscErrorCode ierr = VecDestroy(v);
  checkError(ierr);
}

}
