#include <excafe/util/type_info.hpp>
#include <string>
#include <boost/functional/hash.hpp>

namespace excafe
{

namespace util
{

std::size_t TypeInfo::hashValue() const
{
  return boost::hash<const std::string>()(info->name());
}

std::string TypeInfo::name() const
{
  return info->name();
}

std::size_t hash_value(const TypeInfo& i)
{
  return i.hashValue();
}

}

}
