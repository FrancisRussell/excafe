#ifndef SIMPLE_CFD_UTILITY_HPP
#define SIMPLE_CFD_UTILITY_HPP

#include<iostream>

namespace cfd
{

template<typename T>
class unordered_pair_compare : public std::binary_function<std::pair<T,T>, std::pair<T,T>, bool>
{
public:
  bool operator()(const std::pair<T,T>& p1, const std::pair<T,T>& p2) const
  {
    const std::pair<T,T> sorted_p1 = p1.first < p1.second ? p1 : std::make_pair(p1.second, p1.first);
    const std::pair<T,T> sorted_p2 = p2.first < p2.second ? p2 : std::make_pair(p2.second, p2.first);
    return sorted_p1 < sorted_p2;
  }
};

}


#endif
