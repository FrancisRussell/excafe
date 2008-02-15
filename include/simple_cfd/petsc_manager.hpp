#ifndef SIMPLE_CFD_NUMERIC_PETSC_MANAGER_HPP
#define SIMPLE_CFD_NUMERIC_PETSC_MANAGER_HPP

#include "petsc.h"

namespace cfd
{

class PETScManager
{
private:
  static PETScManager* manager;
  bool initialised;

  PETScManager(const PETScManager& m);
  PETScManager& operator=(const PETScManager& m);
  PETScManager();

public:
  static PETScManager& instance();
  void init(int& argc, char**& argv);
  bool isInitialised() const;
  ~PETScManager();
};

}

#endif
