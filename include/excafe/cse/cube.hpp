#ifndef EXCAFE_CUBE_HPP
#define EXCAFE_CUBE_HPP

#include <ostream>
#include <map>
#include <utility>
#include <algorithm>
#include <boost/operators.hpp>
#include <boost/interprocess/containers/flat_map.hpp>
#include <excafe/util/lazy_copy.hpp>
#include "cse_fwd.hpp"
#include "literal_info.hpp"

namespace excafe
{

namespace cse
{

namespace detail
{
  struct LiteralWriter
  {
    void operator()(std::ostream& o, const int literal) const
    {
      o << "literal(" << literal << ")";
    }
  };
}

class Cube : boost::totally_ordered<Cube,
             boost::additive<Cube
             > >
{
private:
  typedef boost::container::flat_map<unsigned, int> exponent_map_t;
  util::LazyCopy<exponent_map_t> literalExponents; 

  Cube& merge(const Cube& c, bool negate);
  static int minExponent(int a, int b);

public:
  typedef exponent_map_t::value_type     value_type;
  typedef exponent_map_t::iterator       iterator;
  typedef exponent_map_t::const_iterator const_iterator;

  Cube()
  {
  }

  Cube(const unsigned literal)
  {
    literalExponents->insert(std::make_pair(literal, 1));
  }

  Cube(const unsigned literal, const int exponent)
  {
    literalExponents->insert(std::make_pair(literal, exponent));
  }

  template<typename InputIterator>
  Cube(const InputIterator begin, const InputIterator end) : literalExponents(exponent_map_t(begin, end))
  {
  } 

  const_iterator begin() const
  {
    return literalExponents->begin();
  }

  const_iterator end() const
  {
    return literalExponents->end();
  }

  iterator begin()
  {
    return literalExponents->begin();
  }

  iterator end()
  {
    return literalExponents->end();
  }

  std::size_t size() const
  {
    return literalExponents->size();
  }

  bool operator<(const Cube& c) const
  {
    return *literalExponents < *c.literalExponents;
  }

  bool operator==(const Cube& c) const
  {
    return *literalExponents == *c.literalExponents;
  }


  bool contains(const Cube& c) const;
  Cube operator-() const;
  Cube& operator+=(const Cube& c);
  Cube& operator-=(const Cube& c);
  Cube& operator&=(const Cube& c);
  Cube& operator*=(const int exponent);
  void incrementUseCounts(std::map<LiteralInfo, std::size_t>& freqs) const;

  bool isUnit(const NewLiteralCreator& creator) const;
  bool isNumeric(const NewLiteralCreator& creator) const;
  bool hasCoefficient(const NewLiteralCreator& creator) const;
  std::size_t numMultiplies(const NewLiteralCreator& creator) const;

  template<typename literal_writer>
  void write(std::ostream& o, const literal_writer& writer) const
  {
    if (literalExponents->empty())
    {
      o << "1.0";
    }
    else
    {
      for(Cube::const_iterator iter = begin(); iter!=end(); ++iter)
      {
        if (iter != begin())
          o << "*";

        writer(o, iter->first); 
        
        if (iter->second != 1.0)
          o << "^" << iter->second;
      }
    }
  }
};

std::ostream& operator<<(std::ostream& o, const Cube& c);
Cube merge(const Cube& a, const Cube& b, const Cube& c);

}

}

#endif
