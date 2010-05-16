#ifndef SIMPLE_CFD_NUMERIC_VECTOR_HPP
#define SIMPLE_CFD_NUMERIC_VECTOR_HPP

#include <cstddef>
#include <iosfwd>
#include <boost/operators.hpp>
#include "petsc.h"
#include "petscvec.h"

namespace cfd
{

class PETScVector : boost::additive<PETScVector,
                    boost::multipliable<PETScVector, double
                    > >
{
private:
  Vec v;

  void checkError(const PetscErrorCode ierr) const;

public:
  PETScVector(const PETScVector& orig);
  PETScVector(const unsigned rows);
  PETScVector& operator=(const PETScVector& p);
  PETScVector& operator=(const double s);
  PETScVector& operator*=(const double s);
  PETScVector& operator+=(const PETScVector& p);
  PETScVector& operator-=(const PETScVector& p);
  double two_norm() const;
  std::size_t numRows() const;
  void addValues(const unsigned numValues, const int* indices, const double* values);
  void setValues(const unsigned numValues, const int* indices, const double* values);
  void getValues(const unsigned numValues, const int* indices, double* values) const;
  void zero();
  void extractSubvector(PETScVector& dest, const unsigned numValues, const int* indices) const;
  void addSubvector(const PETScVector& source, const unsigned numValues, const int* indices);
  void reciprocal();
  void assemble();
  void print(std::ostream& out) const;
  Vec getPETScHandle() const;
  ~PETScVector();
};

}

#endif
