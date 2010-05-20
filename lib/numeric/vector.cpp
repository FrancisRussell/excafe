#include "simple_cfd/numeric/vector.hpp"
#include <cassert>
#include <numeric>
#include <vector>
#include <ostream>
#include <iostream>
#include <algorithm>
#include "petsc.h"
#include "petscvec.h"

namespace cfd
{

PETScVector::PETScVector(const PETScVector& orig)
{
  PetscErrorCode ierr = VecDuplicate(orig.v, &v);
  checkError(ierr);
  ierr = VecCopy(orig.v, v);
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

PETScVector& PETScVector::operator=(const double s)
{
  const PetscErrorCode ierr = VecSet(v, s);
  checkError(ierr);
  return *this;
}

PETScVector& PETScVector::operator*=(const double s)
{
  const PetscErrorCode ierr = VecScale(v, s);
  checkError(ierr);
  return *this;
}

PETScVector& PETScVector::operator+=(const PETScVector& p)
{
  const PetscErrorCode ierr = VecAXPY(v, 1.0, p.v);
  checkError(ierr);
  return *this;
}

PETScVector& PETScVector::operator-=(const PETScVector& p)
{
  const PetscErrorCode ierr = VecAXPY(v, -1.0, p.v);
  checkError(ierr);
  return *this;
}

PETScVector PETScVector::operator*(const double s) const
{
  PETScVector result(*this);
  result *= s;
  return result;
}

PETScVector PETScVector::operator+(const PETScVector& p) const
{
  PETScVector result(*this);
  result += p;
  return result;
}

PETScVector PETScVector::operator-(const PETScVector& p) const
{
  PETScVector result(*this);
  result -= p;
  return result;
}

double PETScVector::two_norm() const
{
  PetscReal r;
  const PetscErrorCode ierr = VecNorm(v, NORM_2, &r);
  checkError(ierr);
  return r;
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

void PETScVector::zero()
{
  const PetscErrorCode ierr = VecZeroEntries(v);
  checkError(ierr);
}

void PETScVector::reciprocal()
{
  const PetscErrorCode ierr = VecReciprocal(v);
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

void PETScVector::addSubvector(const PETScVector& source, const unsigned numValues, const int* indices)
{
  PetscScalar* data;
  VecGetArray(source.getPETScHandle(), &data);
  addValues(numValues, indices, data);
  VecRestoreArray(source.getPETScHandle(), &data);
}

void PETScVector::print(std::ostream& out = std::cout) const
{
  const size_t rows = numRows();
  std::vector<int> indices(rows);
  std::vector<double> values(rows);

  for(size_t index = 0; index < rows; ++index)
    indices[index] = index;

  getValues(rows, &indices[0], &values[0]);

  out << "(";
  for(size_t index = 0; index < rows; ++index)
    out << values[index] << (index < rows-1 ? ", " : ")");
}

Vec PETScVector::getPETScHandle() const
{
  return v;
}

PETScVector::~PETScVector()
{
  const PetscErrorCode ierr = VecDestroy(v);
  checkError(ierr);
}

}
