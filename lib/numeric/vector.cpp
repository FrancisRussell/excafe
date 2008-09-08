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

std::size_t PETScVector::numRows() const
{
  PetscInt size;
  const PetscErrorCode ierr = VecGetSize(v, &size);
  checkError(ierr);
  return size;
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

void PETScVector::extractSubvector(PETScVector& dest, const unsigned numValues, const int* indices) const
{
  PetscScalar* data;
  VecGetArray(dest.getPETScHandle(), &data);
  getValues(numValues, indices, data);
  VecRestoreArray(dest.getPETScHandle(), &data);
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
