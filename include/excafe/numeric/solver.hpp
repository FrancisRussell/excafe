#ifndef EXCAFE_NUMERIC_SOLVER_HPP
#define EXCAFE_NUMERIC_SOLVER_HPP

#include <string>
#include "excafe_fwd.hpp"
#include "petsc.h"
#include "petscksp.h"
#include "petscpc.h"

namespace excafe
{

class PETScKrylovSolver
{
private:
  KSP ksp;
  PC pc;
  PetscReal rtol;
  PetscReal atol;
  PetscInt maxIts;
  bool preconditionerEnabled;

  void checkError(const PetscErrorCode ierr) const;
  void updateTolerances();
  void updatePreconditioner();

public:
  PETScKrylovSolver();
  void setMaxIterations(const std::size_t maxIter);
  void setRelativeTolerance(const double t);
  void setAbsoluteTolerance(const double t);
  void enablePreconditioner(const bool enable);
  void solve(const PETScMatrix& a, PETScVector& x, const PETScVector& b);
  bool converged() const;
  std::string getConvergedReason() const;
  ~PETScKrylovSolver();
};

}
#endif
