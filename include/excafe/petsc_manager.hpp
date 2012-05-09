#ifndef EXCAFE_NUMERIC_PETSC_MANAGER_HPP
#define EXCAFE_NUMERIC_PETSC_MANAGER_HPP

#include "excafe_fwd.hpp"
#include "petsc.h"

namespace excafe
{

class PETScManager
{
private:
  friend class util::Singleton<PETScManager>;
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
