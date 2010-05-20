#include <simple_cfd/util/singleton.hpp>
#include <simple_cfd/petsc_manager.hpp>
#include <simple_cfd/exception.hpp>

namespace cfd
{

PETScManager::PETScManager() : initialised(false)
{
}

PETScManager& PETScManager::instance()
{
  return util::Singleton<PETScManager>::getInstance();
}

void PETScManager::init(int& argc, char**& argv)
{
  if (!initialised)
  {
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    initialised = true;
  }
  else
  {
    CFD_EXCEPTION("Attempted to initialise PETScManager multiple times.");
  }
}

bool PETScManager::isInitialised() const
{
  return initialised;
}
  
PETScManager::~PETScManager()
{
  PetscFinalize();
}

}

