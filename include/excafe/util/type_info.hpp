#ifndef EXCAFE_UTIL_TYPE_INFO_HPP
#define EXCAFE_UTIL_TYPE_INFO_HPP

#include <typeinfo>
#include <algorithm>
#include <string>
#include <typeinfo>

namespace excafe
{

namespace util
{

class TypeInfo
{
private:
  const std::type_info* info;

public:
  TypeInfo(const std::type_info& _info) : 
    info(&_info)
  {
  }

  bool operator==(const TypeInfo& i) const
  {
    return *info == *i.info;
  }

  bool operator<(const TypeInfo& i) const
  {
    return info->before(*i.info);
  }

  void swap(TypeInfo& i)
  {
    std::swap(info, i.info);
  }

  std::size_t hashValue() const;

  std::string name() const;
};

std::size_t hash_value(const TypeInfo& i);

}

}

#endif
