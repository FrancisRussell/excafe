#ifndef SIMPLE_CFD_NUMERIC_SOLVER_HPP
#define SIMPLE_CFD_NUMERIC_SOLVER_HPP

#include <string>
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
  PetscReal rtol;
  PetscInt maxIts;

  void checkError(const PetscErrorCode ierr) const;
  void updateTolerances();

public:
  PETScKrylovSolver();
  void setMaxIterations(const std::size_t maxIter);
  void setRelativeTolerance(const double t);
  void solve(PETScMatrix& a, PETScVector& x, PETScVector& b);
  bool converged() const;
  std::string getConvergedReason() const;
  ~PETScKrylovSolver();
};

}
#endif
