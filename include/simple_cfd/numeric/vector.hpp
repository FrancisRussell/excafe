#ifndef SIMPLE_CFD_NUMERIC_VECTOR_HPP
#define SIMPLE_CFD_NUMERIC_VECTOR_HPP

#include <cstddef>
#include "petsc.h"
#include "petscvec.h"

namespace cfd
{

class PETScVector
{
private:
  Vec v;

  void checkError(const PetscErrorCode ierr) const;

public:
  PETScVector(const PETScVector& orig);
  PETScVector(const unsigned rows);
  PETScVector& operator=(const PETScVector& p);
  std::size_t numRows() const;
  void addValues(const unsigned numValues, const int* indices, const double* values);
  void setValues(const unsigned numValues, const int* indices, const double* values);
  void getValues(const unsigned numValues, const int* indices, double* values) const;
  void extractSubvector(PETScVector& dest, const unsigned numValues, const int* indices) const;
  void assemble();
  Vec getPETScHandle();
  ~PETScVector();
};

}

#endif
