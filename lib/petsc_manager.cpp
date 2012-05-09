#include <excafe/util/singleton.hpp>
#include <excafe/petsc_manager.hpp>
#include <excafe/exception.hpp>

namespace excafe
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

