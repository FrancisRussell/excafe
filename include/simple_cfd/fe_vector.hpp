#ifndef SIMPLE_CFD_FE_VECTOR_HPP
#define SIMPLE_CFD_FE_VECTOR_HPP

#include <vector>
#include "numeric/vector.hpp"

namespace cfd
{

template<typename C>
class FEVector
{
private:
  typedef C cell_type;
  typedef finite_element<cell_type> finite_element_t;
  typedef typename dof_map<cell_type>::dof_t dof_t;
  const dof_map<cell_type> rowMappings;
  PETScVector vector;
  
  void addOrSetValues(const unsigned rows, const dof_t* rowDofs, const double* values, const bool add)
  {
    std::vector<int> rowIndices(rows);

    for(unsigned row=0; row<rows; ++row)
      rowIndices[row] = rowMappings.getGlobalIndex(rowDofs[row]);

    if (add)
      vector.addValues(rows, &rowIndices[0], values);
    else
      vector.setValues(rows, &rowIndices[0], values);

  }

public:
  FEVector(const dof_map<cell_type>& _rowMappings) : rowMappings(_rowMappings), vector(rowMappings.getDegreesOfFreedomCount())
  {
  }

  FEVector(const dof_map<cell_type>& _rowMappings, const PETScVector& v) : rowMappings(_rowMappings), 
                                                                           vector(v)
  {
  }

  void addValues(const unsigned rows, const dof_t* rowDofs, const double* values)
  {
    addOrSetValues(rows, rowDofs, values, true);
  }

  void setValues(const unsigned rows, const dof_t* rowDofs, const double* values)
  {
    addOrSetValues(rows, rowDofs, values, false);
  }

  void getValues(const unsigned rows, const dof_t* rowDofs, double* values) const
  {
    std::vector<int> rowIndices(rows);

    for(unsigned row=0; row<rows; ++row)
      rowIndices[row] = rowMappings.getGlobalIndex(rowDofs[row]);

    vector.getValues(rows, &rowIndices[0], values);
  }

  void assemble()
  {
    vector.assemble();
  }

  FEVector extractSubvector(const finite_element_t* rowElement) const
  {
    const dof_map<cell_type> newRowMappings = rowMappings.extractDofs(rowElement);
    std::vector<int> rowIndices = newRowMappings.getIndices(rowMappings);
    
    PETScVector subVector = vector.extractSubvector(rowIndices.size(), &rowIndices[0]);
    return FEVector(newRowMappings, subVector);
  }

  PETScVector& getVectorHandle()
  {
    return vector;
  }
};

}

#endif
