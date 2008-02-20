#ifndef SIMPLE_CFD_NUMERIC_SOLVER_HPP
#define SIMPLE_CFD_NUMERIC_SOLVER_HPP

#include "petsc.h"
#include "petscksp.h"
#include "petscpc.h"
#include <cassert>

namespace cfd
{

class PETScKrylovSolver
{
private:
  KSP ksp;
  PC pc;

  void checkError(const PetscErrorCode ierr) const
  {
    assert(ierr == 0);
  }

public:
  PETScKrylovSolver()
  {
    PetscErrorCode ierr;

    ierr = KSPCreate(PETSC_COMM_SELF, &ksp);
    checkError(ierr);

    ierr = KSPGetPC(ksp, &pc);
    checkError(ierr);

    ierr = KSPSetType(ksp, KSPGMRES);
    checkError(ierr);

    ierr = KSPMonitorSet(ksp, KSPMonitorDefault, PETSC_NULL, PETSC_NULL);
    checkError(ierr);

    ierr = PCSetType(pc, PCNONE);
    checkError(ierr);

    //ierr = PCFactorSetAllowDiagonalFill(pc);
    //checkError(ierr);
  }

  void solve(PETScMatrix& a, PETScVector& x, PETScVector& b)
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

  ~PETScKrylovSolver()
  {
    KSPDestroy(ksp);
  }
};

}
#endif
