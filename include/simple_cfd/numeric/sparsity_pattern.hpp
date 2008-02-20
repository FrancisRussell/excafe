#ifndef SIMPLE_CFD_NUMERIC_SPARSITY_PATTERN_HPP
#define SIMPLE_CFD_NUMERIC_SPARSITY_PATTERN_HPP

#include<vector>
#include<set>

namespace cfd
{

class SparsityPattern
{
private:
  const unsigned rows;
  const unsigned cols;
  std::vector< std::set<unsigned> > pattern;

public:
  SparsityPattern(const unsigned r, const unsigned c);
  void insert(const unsigned row, const unsigned col);
  unsigned getRows() const;
  unsigned getCols() const;
  std::vector<int> getNonZerosPerRow() const;
};

}

#endif
