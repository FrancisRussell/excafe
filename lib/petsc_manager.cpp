#include "petsc_manager.hpp"

namespace cfd
{

PETScManager* PETScManager::manager;
  
PETScManager::PETScManager() : initialised(false)
{
}

PETScManager& PETScManager::instance()
{
  if (manager == NULL)
    manager = new PETScManager();

  return *manager;
}

void PETScManager::init(int& argc, char**& argv)
{
  if (!initialised)
  {
    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);
    initialised = true;
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

