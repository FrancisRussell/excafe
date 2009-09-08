#ifndef SIMPLE_CFD_UTIL_TAG_HPP
#define SIMPLE_CFD_UTIL_TAG_HPP

namespace cfd
{

namespace util
{

template<typename T>
class tag
{
public:
  bool operator==(const T& t) const
  {
    return true;
  }

  bool operator<(const T& t) const
  {
    return false;
  }
};

}

}

#endif
