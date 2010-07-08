#include <string>
#include <boost/python/wrapper.hpp>
#include <boost/python/import.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/object.hpp>
#include <boost/python/exec.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <simple_cfd/exception.hpp>
#include <simple_cfd/util/singleton.hpp>
#include <simple_cfd/python_manager.hpp>
#include <simple_cfd/numeric/sympy_bridge.hpp>

namespace cfd
{

void PythonManager::init()
{
  using boost::python::exec;
  using boost::python::import;
  using boost::python::object;
  using boost::python::error_already_set;

  if (!initialised)
  {
    try
    {
      Py_Initialize();

      object main = import("__main__");
      global = main.attr("__dict__");

      // Call the SymPy bridge init code
      cfd::detail::sympy_bridge::init(global);

      initialised = true;
    }
    catch(error_already_set&)
    {
      PyErr_Print();
    }
  }
}

PythonManager::PythonManager() : initialised(false)
{
}

PythonManager& PythonManager::instance()
{
  return util::Singleton<PythonManager>::getInstance();
}

boost::python::object PythonManager::execute(const std::string& code, boost::python::dict& localScope)
{
  if (!initialised)
    CFD_EXCEPTION("PythonManager must be initialised before executing code.");

  return boost::python::exec(code.c_str(), global, localScope);
}

PythonManager::~PythonManager()
{
}

}
