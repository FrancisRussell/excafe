#ifndef EXCAFE_NUMERIC_SPARSITY_PATTERN_HPP
#define EXCAFE_NUMERIC_SPARSITY_PATTERN_HPP

#include<vector>
#include<set>

namespace excafe
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
