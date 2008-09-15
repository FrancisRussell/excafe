#include "simple_cfd/numeric/solver.hpp"
#include "simple_cfd/numeric/matrix.hpp"
#include "simple_cfd/numeric/vector.hpp"
#include "simple_cfd/numeric/sparsity_pattern.hpp"
#include "petsc.h"
#include "petscksp.h"
#include "petscpc.h"
#include <cassert>

namespace cfd
{

void PETScKrylovSolver::checkError(const PetscErrorCode ierr) const
{
  assert(ierr == 0);
}

PETScKrylovSolver::PETScKrylovSolver()
{
  PetscErrorCode ierr;

  ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
  checkError(ierr);
  
  ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  checkError(ierr);

  ierr = KSPGetPC(ksp, &pc);
  checkError(ierr);

  ierr = KSPSetType(ksp, KSPGMRES);
  checkError(ierr);

  ierr = KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
  checkError(ierr);

  ierr = PCSetType(pc, PCNONE);
  checkError(ierr);
}

void PETScKrylovSolver::solve(PETScMatrix& a, PETScVector& x, PETScVector& b)
{
  PetscErrorCode ierr;

  ierr = KSPSetOperators(ksp, a.getPETScHandle(), a.getPETScHandle(), SAME_NONZERO_PATTERN);
  checkError(ierr);

  ierr = KSPSolve(ksp, b.getPETScHandle(), x.getPETScHandle());
  checkError(ierr);

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
    
  if (reason < 0)
    assert(0 && "Solver failed to converge!");
}

PETScKrylovSolver::~PETScKrylovSolver()
{
  KSPDestroy(ksp);
}

}
