#ifndef SMALL_MAP_HPP
#define SMALL_MAP_HPP

#include <functional>
#include <vector>
#include <algorithm>
#include <utility>

namespace cfd
{

namespace util
{

template<typename Key, typename Data, 
  typename Compare = std::less<Key>, 
  typename Alloc = std::allocator< std::pair<Key, Data> > >
class SmallMap
{
public:
  typedef Key key_type;
  typedef Data data_type;
  typedef std::pair<key_type, data_type> value_type;
  typedef Compare key_compare;
  typedef Alloc allocator_type;
 
  class value_compare : public std::binary_function<value_type, value_type, bool>
  {
  private:
    key_compare compare;

  public:
    value_compare(const key_compare& _compare) : compare(_compare) {}
    bool operator()(const value_type& a, const value_type& b) const
    {
      return compare(a.first, b.first);
    }
  };

  typedef value_type* pointer;
  typedef value_type& reference;
  typedef const reference const_reference;
  typedef typename std::vector<value_type>::size_type size_type;
  typedef typename std::vector<value_type>::difference_type difference_type;
  typedef typename std::vector<value_type>::iterator iterator;
  typedef typename std::vector<value_type>::const_iterator const_iterator;
  typedef typename std::vector<value_type>::reverse_iterator reverse_iterator;
  typedef typename std::vector<value_type>::const_reverse_iterator const_reverse_iterator;

private:
  key_compare keyCompare;
  value_compare valueCompare;
  std::vector<value_type, allocator_type> values;

  bool equal(const key_type& a, const key_type& b) const
  {
    return !keyCompare(a,b) && !keyCompare(b, a);
  }

  bool equal(const value_type& a, const value_type& b) const
  {
    return !valueCompare(a,b) && !valueCompare(b, a);
  }

public:
  iterator begin()
  {
    return values.begin();
  }

  iterator end()
  {
    return values.end();
  }

  const_iterator begin() const
  {
    return values.begin();
  }

  const_iterator end() const
  {
    return values.end();
  }

  reverse_iterator rbegin()
  {
    return values.rbegin();
  }

  reverse_iterator rend()
  {
    return values.rend();
  }

  const_reverse_iterator rbegin() const
  {
    return values.rbegin();
  }

  const_reverse_iterator rend() const
  {
    return values.rend();
  }

  size_type size() const
  {
    return values.size();
  }

  size_type max_size() const
  {
    return values.max_size();
  }

  bool empty() const
  {
    return values.empty();
  }

  key_compare key_comp() const
  {
    return keyCompare;
  }

  value_compare value_comp() const
  {
    return valueCompare;
  }

  SmallMap() : valueCompare(keyCompare)
  {
  }

  SmallMap(const key_compare& _keyCompare) : keyCompare(_keyCompare), valueCompare(keyCompare)
  {
  }

  template<typename InputIterator>
  SmallMap(const InputIterator f, const InputIterator l) : valueCompare(keyCompare)
  {
    insert(f, l);
  }

  template<typename InputIterator>
  SmallMap(const InputIterator f, const InputIterator l, const key_compare& _keyCompare) :
    keyCompare(_keyCompare), valueCompare(keyCompare)
  {
    insert(f, l);
  }

  SmallMap(const SmallMap& _map) : keyCompare(_map.keyCompare), valueCompare(_map.valueCompare),
    values(_map.values)
  {
  }

  SmallMap& operator=(const SmallMap& _map)
  {
    keyCompare = _map.keyCompare;
    valueCompare = _map.valueCompare;
    values = _map.values;
    return *this;
  }

  void swap(SmallMap& _map)
  {
    std::swap(keyCompare, _map.keyCompare);
    std::swap(valueCompare, _map.valueCompare);
    std::swap(values, _map.values);
  }

  std::pair<iterator, bool> insert(const value_type& x)
  {
    const iterator pos = lower_bound(x.first);
    if (pos!=end() && equal(x, *pos))
      return std::make_pair(pos, false);
    else
      return std::make_pair(values.insert(pos, x), true);
  }

  iterator insert(const iterator pos, const value_type& x)
  {
    if ((pos==begin() || keyCompare((pos-1)->first, x.first)) && (pos==end() || keyCompare(x.first, pos->first)))
        return values.insert(pos, x);
    
    return insert(x).first;
  }

  template<typename InputIterator>
  void insert(const InputIterator f, const InputIterator l)
  {
    iterator pos = end();
    for(InputIterator current = f; current!=l; ++current)
    {
      pos = insert(pos, *current);
    }
  }

  void erase(const iterator pos)
  {
    values.erase(pos);
  }

  void erase(const iterator first, const iterator last)
  {
    values.erase(first, last);
  }

  void clear()
  {
    values.clear();
  }

  iterator find(const key_type& k)
  {
    const iterator lower = lower_bound(k);
    if (lower == end() || !equal(lower->first, k))
      return end();
    else
      return lower;
  }

  const_iterator find(const key_type& k) const
  {
    const const_iterator lower = lower_bound(k);
    if (lower == end() || !equal(lower->first, k))
      return end();
    else
      return lower;
  }

  size_type count(const key_type& k) const
  {
    const std::pair<const_iterator, const_iterator> range = equal_range(k);
    return range.second - range.first;
  }

  iterator lower_bound(const key_type& k)
  {
    return std::lower_bound(begin(), end(), value_type(k, data_type()), valueCompare);
  }

  const_iterator lower_bound(const key_type& k) const
  {
    return std::lower_bound(begin(), end(), value_type(k, data_type()), valueCompare);
  }

  iterator upper_bound(const key_type& k)
  {
    return std::upper_bound(begin(), end(), value_type(k, data_type()), valueCompare);
  }

  const_iterator upper_bound(const key_type& k) const
  {
    return std::upper_bound(begin(), end(), value_type(k, data_type()), valueCompare);
  }

  std::pair<iterator, iterator> equal_range(const key_type& k)
  {
    return std::equal_range(begin(), end(), value_type(k, data_type()), valueCompare);
  }

  std::pair<const_iterator, const_iterator> equal_range(const key_type& k) const
  {
    return std::equal_range(begin(), end(), value_type(k, data_type()), valueCompare);
  }

  data_type& operator[](const key_type& k)
  {
    const iterator lower = lower_bound(k);
    if (lower == end() || !equal(lower->first, k))
      return insert(lower, data_type())->second;
    else
      return lower->second;
  }

  bool operator==(const SmallMap& _map) const
  {
    return values == _map.values;
  }

  bool operator<(const SmallMap& _map) const
  {
    return values < _map.values;
  }
};

}

}

namespace std
{

template<typename K, typename D, typename C, typename A>
void swap(cfd::util::SmallMap<K,D,C,A>& a, cfd::util::SmallMap<K,D,C,A>& b)
{
  a.swap(b);
}

}

#endif
