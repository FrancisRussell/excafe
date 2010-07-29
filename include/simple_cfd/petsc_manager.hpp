#ifndef SIMPLE_CFD_NUMERIC_PETSC_MANAGER_HPP
#define SIMPLE_CFD_NUMERIC_PETSC_MANAGER_HPP

#include "simple_cfd_fwd.hpp"
#include "petsc.h"

namespace cfd
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
