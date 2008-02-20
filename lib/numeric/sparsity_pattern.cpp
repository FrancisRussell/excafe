#include "simple_cfd/numeric/sparsity_pattern.hpp"

namespace cfd
{

SparsityPattern::SparsityPattern(const unsigned r, const unsigned c) : rows(r), cols(c), pattern(rows)
{
}

void SparsityPattern::insert(const unsigned row, const unsigned col)
{
  pattern[row].insert(col);
}

unsigned SparsityPattern::getRows() const
{
  return rows;
}

unsigned SparsityPattern::getCols() const
{
  return cols;
}

std::vector<int> SparsityPattern::getNonZerosPerRow() const
{
  std::vector<int> nnz;
  nnz.reserve(rows);

  for(std::vector< std::set<unsigned> >::const_iterator rowIter(pattern.begin()); rowIter!=pattern.end(); ++rowIter)
    nnz.push_back(rowIter->size());

  return nnz;
}

}
