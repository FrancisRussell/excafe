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

PETScKrylovSolver::PETScKrylovSolver() : rtol(PETSC_DEFAULT), atol(PETSC_DEFAULT), maxIts(PETSC_DEFAULT), preconditionerEnabled(true)
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

  ierr = PCSetType(pc, PCSOR);
  checkError(ierr);

  updateTolerances();
  updatePreconditioner();
}

void PETScKrylovSolver::solve(const PETScMatrix& a, PETScVector& x, const PETScVector& b)
{
  PetscErrorCode ierr;

  ierr = KSPSetOperators(ksp, a.getPETScHandle(), a.getPETScHandle(), SAME_NONZERO_PATTERN);
  checkError(ierr);

  ierr = KSPSolve(ksp, b.getPETScHandle(), x.getPETScHandle());
  checkError(ierr);
}

bool PETScKrylovSolver::converged() const
{
  KSPConvergedReason reason;
  const PetscErrorCode ierr = KSPGetConvergedReason(ksp, &reason);
  checkError(ierr);
  return reason > 0;
}

std::string PETScKrylovSolver::getConvergedReason() const
{
  KSPConvergedReason reason;
  const PetscErrorCode ierr = KSPGetConvergedReason(ksp, &reason);
  checkError(ierr);

  switch(reason)
  {
    case KSP_CONVERGED_RTOL:            return "KSP_CONVERGED_RTOL";
    case KSP_CONVERGED_ATOL:            return "KSP_CONVERGED_ATOL";
    case KSP_CONVERGED_ITS:             return "KSP_CONVERGED_ITS";
    case KSP_CONVERGED_CG_NEG_CURVE:    return "KSP_CONVERGED_CG_NEG_CURVE";
    case KSP_CONVERGED_CG_CONSTRAINED:  return "KSP_CONVERGED_CG_CONSTRAINED";
    case KSP_CONVERGED_STEP_LENGTH:     return "KSP_CONVERGED_STEP_LENGTH";
    case KSP_CONVERGED_HAPPY_BREAKDOWN: return "KSP_CONVERGED_HAPPY_BREAKDOWN";
    case KSP_DIVERGED_NULL:             return "KSP_DIVERGED_NULL";
    case KSP_DIVERGED_ITS:              return "KSP_DIVERGED_ITS";
    case KSP_DIVERGED_DTOL:             return "KSP_DIVERGED_DTOL";
    case KSP_DIVERGED_BREAKDOWN:        return "KSP_DIVERGED_BREAKDOWN";
    case KSP_DIVERGED_BREAKDOWN_BICG:   return "KSP_DIVERGED_BREAKDOWN_BICG";
    case KSP_DIVERGED_NONSYMMETRIC:     return "KSP_DIVERGED_NONSYMMETRIC";
    case KSP_DIVERGED_INDEFINITE_PC:    return "KSP_DIVERGED_INDEFINITE_PC";
    case KSP_DIVERGED_NAN:              return "KSP_DIVERGED_NAN";
    case KSP_DIVERGED_INDEFINITE_MAT:   return "KSP_DIVERGED_INDEFINITE_MAT";
    case KSP_CONVERGED_ITERATING:       return "KSP_CONVERGED_ITERATING";
    default:                            return "Unknown reason";
  }
}

void PETScKrylovSolver::updateTolerances()
{
  const PetscErrorCode ierr = KSPSetTolerances(ksp, rtol, atol, PETSC_DEFAULT, maxIts);
  checkError(ierr);
}

void PETScKrylovSolver::updatePreconditioner()
{
  PetscErrorCode ierr = PCSetType(pc, preconditionerEnabled ? PCSOR : PCNONE);
  checkError(ierr);
}

void PETScKrylovSolver::enablePreconditioner(const bool enable)
{
  preconditionerEnabled = enable;
  updatePreconditioner();
}

void PETScKrylovSolver::setMaxIterations(const std::size_t maxIter)
{
  maxIts = maxIter;
  updateTolerances();
}

void PETScKrylovSolver::setRelativeTolerance(const double t)
{
  rtol = t;
  updateTolerances();
}

void PETScKrylovSolver::setAbsoluteTolerance(const double t)
{
  atol = t;
  updateTolerances();
}

PETScKrylovSolver::~PETScKrylovSolver()
{
  KSPDestroy(&ksp);
}

}
