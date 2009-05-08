#ifndef SIMPLE_CFD_FE_VECTOR_HPP
#define SIMPLE_CFD_FE_VECTOR_HPP

#include <vector>
#include <ostream>
#include <iostream>
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

  FEVector(const FEVector& v) : rowMappings(v.rowMappings), vector(v.vector)
  {
  }

  dof_map<cell_type> getRowMappings() const
  {
    return rowMappings;
  }

  FEVector& operator=(const FEVector& f)
  {
    assert(rowMappings == f.rowMappings);
    vector = f.vector;
    return *this;
  }

  FEVector& operator=(const double s)
  {
    vector = s;
    return *this;
  }

  FEVector& operator=(const double s)
  {
    vector = s;
  }

  FEVector& operator*=(const double d)
  {
    vector *= d;
    return *this;
  }

  FEVector& operator+=(const FEVector& f)
  {
    assert(rowMappings == f.rowMappings);
    vector += f.vector;
    return *this;
  }
  
  FEVector& operator-=(const FEVector& f)
  {
    assert(rowMappings == f.rowMappings);
    vector -= f.vector;
    return *this;
  }

  FEVector operator*(const double d) const
  {
    return FEVector(rowMappings, vector * d);
  }

  FEVector operator+(const FEVector& f) const
  {
    assert(rowMappings == f.rowMappings);
    return FEVector(rowMappings, vector + f.vector);
  }

  FEVector operator-(const FEVector& f) const
  {
    assert(rowMappings == f.rowMappings);
    return FEVector(rowMappings, vector - f.vector);
  }

  double two_norm() const
  {
    return vector.two_norm();
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

  void zero()
  {
    vector.zero();
  }

  void reciprocal()
  {
    vector.reciprocal();
  }

  void assemble()
  {
    vector.assemble();
  }

  void extractSubvector(FEVector& s) const
  {
    const std::vector<int> rowIndices = s.rowMappings.getIndices(rowMappings);
    vector.extractSubvector(s.vector, rowIndices.size(), &rowIndices[0]);
  }

  void addSubvector(const FEVector& s)
  {
    const std::vector<int> rowIndices = s.rowMappings.getIndices(rowMappings);
    vector.addSubvector(s.vector, rowIndices.size(), &rowIndices[0]);
  }

  void print(std::ostream& out = std::cout)
  {
    vector.print(out);
  }

  PETScVector& getVectorHandle()
  {
    return vector;
  }
};

}

#endif
