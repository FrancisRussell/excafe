#ifndef SIMPLE_CFD_FE_VECTOR_HPP
#define SIMPLE_CFD_FE_VECTOR_HPP

#include <vector>
#include <ostream>
#include <iostream>
#include "numeric/vector.hpp"
#include "dof_map.hpp"

namespace cfd
{

template<std::size_t D>
class DiscreteField
{
private:
  static const std::size_t dimension = D;
  typedef FiniteElement<dimension> finite_element_t;
  typedef typename DofMap<dimension>::dof_t dof_t;
  DofMap<dimension> rowMappings;
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
  DiscreteField(const DofMap<dimension>& _rowMappings) : rowMappings(_rowMappings), vector(rowMappings.getDegreesOfFreedomCount())
  {
  }

  DiscreteField(const DofMap<dimension>& _rowMappings, const PETScVector& v) : rowMappings(_rowMappings), 
                                                                           vector(v)
  {
  }

  DiscreteField(const DiscreteField& v) : rowMappings(v.rowMappings), vector(v.vector)
  {
  }

  void swap(DiscreteField& f)
  {
    rowMappings.swap(f.rowMappings);
    vector.swap(f.vector);
  }

  bool isComposite() const
  {
    return rowMappings.isComposite();
  }

  std::size_t getRank() const
  {
    assert(!isComposite());
    return rowMappings.getFiniteElement()->getRank();
  }

  std::size_t getDimension() const
  {
    assert(!isComposite());
    return rowMappings.getFiniteElement()->getDimension();
  }

  const finite_element_t* getElement() const
  {
    assert(!isComposite());
    return rowMappings.getFiniteElement();
  }

  DofMap<dimension> getRowMappings() const
  {
    return rowMappings;
  }

  DiscreteField& operator=(const DiscreteField& f)
  {
    assert(rowMappings == f.rowMappings);
    vector = f.vector;
    return *this;
  }

  DiscreteField& operator=(const double s)
  {
    vector = s;
    return *this;
  }

  DiscreteField& operator*=(const double d)
  {
    vector *= d;
    return *this;
  }

  DiscreteField& operator+=(const DiscreteField& f)
  {
    assert(rowMappings == f.rowMappings);
    vector += f.vector;
    return *this;
  }
  
  DiscreteField& operator-=(const DiscreteField& f)
  {
    assert(rowMappings == f.rowMappings);
    vector -= f.vector;
    return *this;
  }

  DiscreteField operator*(const double d) const
  {
    return DiscreteField(rowMappings, vector * d);
  }

  DiscreteField operator+(const DiscreteField& f) const
  {
    assert(rowMappings == f.rowMappings);
    return DiscreteField(rowMappings, vector + f.vector);
  }

  DiscreteField operator-(const DiscreteField& f) const
  {
    assert(rowMappings == f.rowMappings);
    return DiscreteField(rowMappings, vector - f.vector);
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

  void extractField(DiscreteField& s) const
  {
    const std::vector<int> rowIndices = s.rowMappings.getIndices(rowMappings);
    vector.extractSubvector(s.vector, rowIndices.size(), &rowIndices[0]);
    s.assemble();
  }

  void addField(const DiscreteField& s)
  {
    const std::vector<int> rowIndices = s.rowMappings.getIndices(rowMappings);
    vector.addSubvector(s.vector, rowIndices.size(), &rowIndices[0]);
    assemble();
  }

  void print(std::ostream& out = std::cout)
  {
    vector.print(out);
  }

  PETScVector& getVectorHandle()
  {
    return vector;
  }

  const PETScVector& getVectorHandle() const
  {
    return vector;
  }

  DiscreteField project(const DofMap<dimension>& newDofMap)
  {
    const DofMap<dimension> intermediateMap(rowMappings.intersect(newDofMap));
    
    DiscreteField result(newDofMap);
    DiscreteField intermediateField(intermediateMap);

    extractField(intermediateField);
    result.addField(intermediateField);

    return result;
  }
};

}

#endif
