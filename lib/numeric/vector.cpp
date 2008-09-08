#include "simple_cfd/numeric/vector.hpp"
#include <cassert>
#include "petsc.h"
#include "petscvec.h"

namespace cfd
{

PETScVector::PETScVector(const PETScVector& orig)
{
  const PetscErrorCode ierr = VecDuplicate(orig.v, &v);
  checkError(ierr);
}

void PETScVector::checkError(const PetscErrorCode ierr) const
{
  assert(ierr == 0);
}

PETScVector::PETScVector(const unsigned rows)
{
  const PetscErrorCode ierr = VecCreateSeq(PETSC_COMM_SELF, rows, &v);
  checkError(ierr);
}

PETScVector& PETScVector::operator=(const PETScVector& p)
{
  const PetscErrorCode ierr = VecCopy(p.v, v);
  checkError(ierr);
  return *this;
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

PETScVector PETScVector::extractSubvector(const unsigned numValues, const int* indices) const
{
  PETScVector result(numValues);
  PetscScalar* data;

  VecGetArray(result.getPETScHandle(), &data);
  getValues(numValues, indices, data);
  VecRestoreArray(result.getPETScHandle(), &data);

  return result;
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
