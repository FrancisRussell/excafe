#ifndef SIMPLE_CFD_NUMERIC_SOLVER_HPP
#define SIMPLE_CFD_NUMERIC_SOLVER_HPP

#include "simple_cfd_fwd.hpp"
#include "petsc.h"
#include "petscksp.h"
#include "petscpc.h"

namespace cfd
{

class PETScKrylovSolver
{
private:
  KSP ksp;
  PC pc;

  void checkError(const PetscErrorCode ierr) const;

public:
  PETScKrylovSolver();
  void solve(PETScMatrix& a, PETScVector& x, PETScVector& b);
  ~PETScKrylovSolver();
};

}
#endif
