#ifndef EXCAFE_UTIL_TAG_HPP
#define EXCAFE_UTIL_TAG_HPP

namespace excafe
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
