#include <string>
#include <boost/python/wrapper.hpp>
#include <boost/python/import.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/object.hpp>
#include <boost/python/exec.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <excafe/exception.hpp>
#include <excafe/util/singleton.hpp>
#include <excafe/python_manager.hpp>
#include <excafe/numeric/sympy_bridge.hpp>

namespace excafe
{

void PythonManager::init()
{
  using boost::python::exec;
  using boost::python::import;
  using boost::python::object;
  using boost::python::error_already_set;

  if (!Py_IsInitialized())
  {
    try
    {
      Py_Initialize();

      object main = import("__main__");
      global = main.attr("__dict__");

      // Call the SymPy bridge init code
      excafe::detail::sympy_bridge::init(global);
    }
    catch(error_already_set&)
    {
      PyErr_Print();
    }
  }
  else
  {
    CFD_EXCEPTION("Double Python initialisation detected.");
  }
}

PythonManager::PythonManager()
{
}

PythonManager& PythonManager::instance()
{
  return util::Singleton<PythonManager>::getInstance();
}

boost::python::object PythonManager::getObject(const std::string& name)
{
  return global[name];
}

boost::python::object PythonManager::execute(const std::string& code, boost::python::dict& localScope)
{
  if (!Py_IsInitialized())
    CFD_EXCEPTION("PythonManager must be initialised before executing code.");

  return boost::python::exec(code.c_str(), global, localScope);
}

PythonManager::~PythonManager()
{
  // Enable this when Boost says it's OK
  // Py_Finalize();
}

}
