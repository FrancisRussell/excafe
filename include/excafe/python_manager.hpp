#ifndef EXCAFE_PYTHON_MANAGER_HPP
#define EXCAFE_PYTHON_MANAGER_HPP

#include <string>
#include <boost/python/object.hpp>
#include <boost/python/dict.hpp>
#include <boost/noncopyable.hpp>
#include <excafe/util/singleton.hpp>

namespace excafe
{

class PythonManager : boost::noncopyable
{
private:
  friend class util::Singleton<PythonManager>;
  PythonManager();

  bool initialised;
  boost::python::object global;

public:
  static PythonManager& instance();
  void init();
  boost::python::object getObject(const std::string& name);
  boost::python::object execute(const std::string& code, boost::python::dict& localScope);
  ~PythonManager();
};

}

#endif
